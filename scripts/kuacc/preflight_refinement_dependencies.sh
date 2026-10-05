#!/usr/bin/env bash
# Bounded KUACC dependency preflight; this never submits a production job.
set -eo pipefail

ROOT=${1:?usage: preflight_refinement_dependencies.sh REMOTE_PREFLIGHT_ROOT}
TOOLS="$ROOT/tools"
FIBER="$TOOLS/fiberdock"
MP="$TOOLS/multiprot.Linux"
LEGACY="${MULTIPROT_LEGACY_ROOT:-$TOOLS/multiprot_legacy_runtime}"
REPORT="$ROOT/manifest"
LOGS="$ROOT/logs"
mkdir -p "$REPORT" "$LOGS"

source /etc/bashrc
set -u
module load rosetta/2022.42

ROSETTA_PREPACK=$(command -v docking_prepack_protocol.static.linuxgccrelease || true)
ROSETTA_DOCK=$(command -v docking_protocol.static.linuxgccrelease || true)
ROSETTA_DB=${PRISM_ROSETTA_DB:-/kuacc/apps/rosetta/rosetta_bin_linux_2022.42_bundle/main/database}

manifest="$REPORT/dependency_manifest.tsv"
printf 'timestamp\tcomponent\tpath\tstatus\tdetail\n' > "$manifest"
record() {
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date -Is)" "$1" "$2" "$3" "$4" >> "$manifest"
}
require_file() {
  local component=$1 path=$2
  if [[ -f "$path" ]]; then
    record "$component" "$path" present "$(stat -c '%s bytes mode=%a' "$path")"
  else
    record "$component" "$path" missing "required file"
    return 1
  fi
}
require_exec() {
  local component=$1 path=$2
  if [[ -x "$path" ]]; then
    record "$component" "$path" executable "$(file -b "$path")"
  else
    record "$component" "$path" not_executable "required executable"
    return 1
  fi
}

echo "host=$(hostname)" | tee "$REPORT/host.txt"
echo "rosetta_module=rosetta/2022.42" | tee "$REPORT/rosetta_module.txt"
module show rosetta/2022.42 > "$REPORT/rosetta_module_show.txt" 2>&1 || true
printf 'ROSETTA_PREPACK=%s\nROSETTA_DOCK=%s\nROSETTA_DB=%s\n' \
  "$ROSETTA_PREPACK" "$ROSETTA_DOCK" "$ROSETTA_DB" | tee "$REPORT/rosetta_paths.txt"

[[ -n "$ROSETTA_PREPACK" && -x "$ROSETTA_PREPACK" ]] || { record rosetta "$ROSETTA_PREPACK" missing "prepack executable"; exit 2; }
[[ -n "$ROSETTA_DOCK" && -x "$ROSETTA_DOCK" ]] || { record rosetta "$ROSETTA_DOCK" missing "docking executable"; exit 2; }
[[ -d "$ROSETTA_DB" ]] || { record rosetta "$ROSETTA_DB" missing "database directory"; exit 2; }
record rosetta "$ROSETTA_PREPACK" executable "$(file -b "$ROSETTA_PREPACK")"
record rosetta "$ROSETTA_DOCK" executable "$(file -b "$ROSETTA_DOCK")"
record rosetta "$ROSETTA_DB" present "database directory"

require_exec multiprot "$MP"
require_file multiprot_legacy "$LEGACY/prism.py"
require_file multiprot_legacy "$LEGACY/run_files/mainController.py"
require_file multiprot_legacy "$LEGACY/run_files/structuralAlignment.py"
require_file multiprot_legacy "$LEGACY/template_default"
require_file multiprot_legacy "$LEGACY/template_list"

require_exec fiberdock "$FIBER/FiberDock"
require_exec fiberdock "$FIBER/nma"
require_exec fiberdock "$FIBER/reduce"
require_file fiberdock "$FIBER/buildFiberDockParams.pl"
require_file fiberdock "$FIBER/addHydrogens.pl"
require_file fiberdock "$FIBER/reduce_het_dict.txt"
require_file fiberdock "$FIBER/zero-trial"
for library in Names.CHARMM.db bbdep02.May.sortlib chem.lib exclusions par_all27_prot_na.prm top_all27_prot_na.rtf; do
  require_file fiberdock "$FIBER/lib/$library"
done

sha256sum "$MP" "$FIBER/FiberDock" "$FIBER/nma" "$FIBER/reduce" \
  "$LEGACY/prism.py" "$LEGACY/run_files/mainController.py" > "$REPORT/tool_hashes.sha256"

set +e
timeout 15s "$MP" --help > "$LOGS/multiprot_help.stdout" 2> "$LOGS/multiprot_help.stderr"
mp_rc=$?
set -e
if grep -q '^Usage:' "$LOGS/multiprot_help.stderr"; then
  record multiprot "$MP" smoke_pass "usage emitted; rc=$mp_rc"
else
  record multiprot "$MP" smoke_fail "usage not emitted; rc=$mp_rc"
  exit 3
fi

set +e
timeout 20s "$ROSETTA_PREPACK" -help > "$LOGS/rosetta_prepack_help.stdout" 2> "$LOGS/rosetta_prepack_help.stderr"
prepack_rc=$?
timeout 20s "$ROSETTA_DOCK" -help > "$LOGS/rosetta_dock_help.stdout" 2> "$LOGS/rosetta_dock_help.stderr"
dock_rc=$?
set -e
if grep -Eiq 'Rosetta|database|options|Usage' "$LOGS/rosetta_prepack_help.stdout" "$LOGS/rosetta_prepack_help.stderr"; then
  record rosetta "$ROSETTA_PREPACK" smoke_pass "help emitted; rc=$prepack_rc"
else
  record rosetta "$ROSETTA_PREPACK" smoke_fail "help not emitted; rc=$prepack_rc"
fi
if grep -Eiq 'Rosetta|database|options|Usage' "$LOGS/rosetta_dock_help.stdout" "$LOGS/rosetta_dock_help.stderr"; then
  record rosetta "$ROSETTA_DOCK" smoke_pass "help emitted; rc=$dock_rc"
else
  record rosetta "$ROSETTA_DOCK" smoke_fail "help not emitted; rc=$dock_rc"
fi

set +e
perl -c "$FIBER/buildFiberDockParams.pl" > "$LOGS/buildFiberDockParams.perl-c" 2>&1
params_rc=$?
perl -c "$FIBER/addHydrogens.pl" > "$LOGS/addHydrogens.perl-c" 2>&1
hydrogen_rc=$?
set -e
[[ "$params_rc" -eq 0 ]] && record fiberdock "$FIBER/buildFiberDockParams.pl" syntax_pass "perl -c" || record fiberdock "$FIBER/buildFiberDockParams.pl" syntax_fail "perl -c rc=$params_rc"
[[ "$hydrogen_rc" -eq 0 ]] && record fiberdock "$FIBER/addHydrogens.pl" syntax_pass "perl -c" || record fiberdock "$FIBER/addHydrogens.pl" syntax_fail "perl -c rc=$hydrogen_rc"

ldd "$FIBER/FiberDock" > "$LOGS/FiberDock.ldd" 2>&1 || true
if grep -q 'not found' "$LOGS/FiberDock.ldd"; then
  record fiberdock "$FIBER/FiberDock" link_fail "missing shared library"
else
  record fiberdock "$FIBER/FiberDock" link_pass "ldd has no missing libraries"
fi
for tool in FiberDock nma reduce; do
  set +e
  (cd "$FIBER" && timeout 15s "./$tool" > "$LOGS/${tool}_usage.stdout" 2> "$LOGS/${tool}_usage.stderr")
  tool_rc=$?
  set -e
  if grep -Eiq 'Usage:|arguments:|input.pdb' "$LOGS/${tool}_usage.stdout" "$LOGS/${tool}_usage.stderr"; then
    record fiberdock "$FIBER/$tool" smoke_pass "usage emitted; rc=$tool_rc"
  else
    record fiberdock "$FIBER/$tool" smoke_fail "usage not emitted; rc=$tool_rc"
  fi
done

printf 'preflight_completed=%s\n' "$(date -Is)" | tee "$REPORT/status.txt"
