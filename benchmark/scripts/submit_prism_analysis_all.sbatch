#!/bin/bash
# Submit with:
#   sbatch benchmark/scripts/submit_prism_analysis_all.sbatch
#
# Optional overrides:
#   sbatch --export=ALL,SETS=rigid,medium,difficult,SCORE_TIMEOUT_SEC=5 benchmark/scripts/submit_prism_analysis_all.sbatch

#SBATCH --job-name=prism-bench
#SBATCH --account=ai
#SBATCH --partition=ai
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=benchmark/prism_processed/results/slurm-%J.log

set -euo pipefail

ROOT_DIR="/scratch/rshadi25/GitHub/PRISM"
PRISM_GDRIVE_FILE_ID="${PRISM_GDRIVE_FILE_ID:-1znDfQGJ1SpKnhAX-zzAUC5M8KVg71qgm}"
PRISM_RAW_ARCHIVE="${PRISM_RAW_ARCHIVE:-$ROOT_DIR/benchmark/prism_processed/prism_raw.zip}"
PRISM_RAW_DIR="${PRISM_RAW_DIR:-$ROOT_DIR/benchmark/prism_processed/prism_raw}"
ROSETTA_ROOT="${ROSETTA_ROOT:-$PRISM_RAW_DIR/rosetta_output_1}"
FIXED_ROOT="${FIXED_ROOT:-$ROOT_DIR/benchmark/prism_processed/prism_raw_chainfixed}"
RESULTS_ROOT="${RESULTS_ROOT:-$ROOT_DIR/benchmark/prism_processed/results}"
DATA_DIR="${DATA_DIR:-$ROOT_DIR/benchmark/data}"
SCORE_PYTHON="${SCORE_PYTHON:-$ROOT_DIR/benchmark/prism_processed/env/prism_score_env/bin/python}"
SETS="${SETS:-rigid,medium,difficult}"
METRICS="${METRICS:-dockq,irmsd}"
SCORE_TIMEOUT_SEC="${SCORE_TIMEOUT_SEC:-5}"
VERBOSE_EVERY="${VERBOSE_EVERY:-50}"

CHAIN_FIX_MAP="$RESULTS_ROOT/chain_fix_map.csv"

mkdir -p "$RESULTS_ROOT" "$PRISM_RAW_DIR"

echo "[info] root: $ROOT_DIR"
echo "[info] prism_gdrive_file_id: $PRISM_GDRIVE_FILE_ID"
echo "[info] prism_raw_archive: $PRISM_RAW_ARCHIVE"
echo "[info] prism_raw_dir: $PRISM_RAW_DIR"
echo "[info] rosetta_root: $ROSETTA_ROOT"
echo "[info] fixed_root: $FIXED_ROOT"
echo "[info] results_root: $RESULTS_ROOT"
echo "[info] score_python: $SCORE_PYTHON"
echo "[info] sets: $SETS"
echo "[info] metrics: $METRICS"
echo "[info] score_timeout_sec: $SCORE_TIMEOUT_SEC"

cd "$ROOT_DIR"

download_prism_archive() {
  if [[ -s "$PRISM_RAW_ARCHIVE" ]]; then
    echo "[info] using existing archive: $PRISM_RAW_ARCHIVE"
    return 0
  fi

  mkdir -p "$(dirname "$PRISM_RAW_ARCHIVE")"

  if command -v gdown >/dev/null 2>&1; then
    echo "[info] downloading archive with gdown"
    gdown --fuzzy "https://drive.google.com/file/d/$PRISM_GDRIVE_FILE_ID/view?usp=sharing" -O "$PRISM_RAW_ARCHIVE"
    return 0
  fi

  if command -v curl >/dev/null 2>&1; then
    echo "[info] downloading archive with curl direct link fallback"
    curl -L "https://drive.google.com/uc?export=download&id=$PRISM_GDRIVE_FILE_ID" -o "$PRISM_RAW_ARCHIVE"
    return 0
  fi

  if command -v wget >/dev/null 2>&1; then
    echo "[info] downloading archive with wget direct link fallback"
    wget --no-check-certificate -O "$PRISM_RAW_ARCHIVE" \
      "https://drive.google.com/uc?export=download&id=$PRISM_GDRIVE_FILE_ID"
    return 0
  fi

  echo "[error] no downloader available; install gdown, curl, or wget" >&2
  exit 1
}

extract_prism_archive() {
  if [[ -d "$ROSETTA_ROOT" ]]; then
    echo "[info] using existing rosetta root: $ROSETTA_ROOT"
    return 0
  fi

  rm -rf "$PRISM_RAW_DIR"
  mkdir -p "$PRISM_RAW_DIR"

  case "$PRISM_RAW_ARCHIVE" in
    *.zip)
      unzip -q "$PRISM_RAW_ARCHIVE" -d "$PRISM_RAW_DIR"
      ;;
    *.tar.gz|*.tgz)
      tar -xzf "$PRISM_RAW_ARCHIVE" -C "$PRISM_RAW_DIR"
      ;;
    *.tar)
      tar -xf "$PRISM_RAW_ARCHIVE" -C "$PRISM_RAW_DIR"
      ;;
    *)
      echo "[error] unsupported archive format: $PRISM_RAW_ARCHIVE" >&2
      exit 1
      ;;
  esac

  if [[ -d "$PRISM_RAW_DIR/rosetta_output_1" ]]; then
    ROSETTA_ROOT="$PRISM_RAW_DIR/rosetta_output_1"
    return 0
  fi

  local detected_root
  detected_root="$(find "$PRISM_RAW_DIR" -type d -name rosetta_output_1 | head -n 1 || true)"
  if [[ -n "$detected_root" ]]; then
    ROSETTA_ROOT="$detected_root"
    echo "[info] detected rosetta root: $ROSETTA_ROOT"
    return 0
  fi

  echo "[error] could not find rosetta_output_1 after extracting $PRISM_RAW_ARCHIVE" >&2
  exit 1
}

download_prism_archive
extract_prism_archive

echo "[info] resolved rosetta_root: $ROSETTA_ROOT"

python3 -u benchmark/scripts/fix_model_chain_names.py \
  --rosetta-root "$ROSETTA_ROOT" \
  --fixed-root "$FIXED_ROOT" \
  --map-csv "$CHAIN_FIX_MAP"

python3 -u benchmark/scripts/analyze_prism_all_benchmarks.py \
  --rosetta-root "$ROSETTA_ROOT" \
  --data-dir "$DATA_DIR" \
  --results-root "$RESULTS_ROOT" \
  --score-python "$SCORE_PYTHON" \
  --metrics "$METRICS" \
  --sets "$SETS" \
  --model-fix-map-csv "$CHAIN_FIX_MAP" \
  --score-timeout-sec "$SCORE_TIMEOUT_SEC" \
  --verbose-every "$VERBOSE_EVERY" \
  --download-native-missing

echo "[done] benchmark analysis finished"
