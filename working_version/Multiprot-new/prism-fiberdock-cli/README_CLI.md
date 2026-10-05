# PRISM — FiberDock/MultiProt CLI build

Command-line version of the **original FiberDock + MultiProt** PRISM pipeline
(the algorithm behind the live web server), with the PDB-fetch layer modernised
from the newer TMalign/Rosetta build (`prism-fixed-vFatma`).

## What was changed vs. the web-server source (`prism-main`)

Only the fetch/IO layer and the web-server plumbing were touched. The scientific
pipeline (surface extraction → **MultiProt** alignment → transformation filtering →
**FiberDock** flexible refinement) is unchanged.

| File | Change |
|------|--------|
| `run_files/pdbDownload.py` | FTP (`ftp.wwpdb.org`) → **HTTPS** (`files.pdbj.org`); added **mmCIF fallback** with BeEM→PDB conversion for entries that have no legacy PDB file. |
| `run_files/templateGenerator.py` | Same HTTPS + mmCIF/BeEM fetch grafted into its `pdbDownload()`; template-generation logic untouched. Added the imports the original was missing. |
| `run_files/merge_bundles.py` | **New** helper (from vFatma) — merges BeEM multi-bundle PDBs into one file. |
| `external_tools/BeEM-master/` | **New** tool — mmCIF→PDB converter (compile on target host). |
| `run_files/mainController.py` | MySQL, mail, and HTML-progress calls disabled; keeps MultiProt aligner + FiberDock refiner. Per-job cap lifted 10 → 100. Intermediate folders kept (not deleted). |
| `run_files/structuralAlignment.py` | Removed the MySQL `multiprot` cache; every alignment is now written to a per-job pickle file (`alignment/<interface>_<chain>_<target>`). |
| `run_files/transformationFiltering.py` | Reads those alignment files instead of the MySQL cache. |

`mysqlWriter.py`, `databaseChecker.py`, `sendMail.py`, `startMail.py`,
`htmlWriter.py` are left in the tree but are **no longer imported** by the CLI path.

## One-time setup (on the KU Linux cluster)

Requires **Python 2.7**. No MySQL, no mail server needed.

```bash
cd external_tools

# unpack the FiberDock/MultiProt toolchain
unzip fiberdock.zip
unzip multiprot.zip
unzip naccess.zip
unzip pops.zip          # only if you set external_tool_choice = 0

# naccess needs compiling
cd naccess && csh install.scr && cd ..

# BeEM (mmCIF -> PDB converter)
cd BeEM-master && g++ -O3 BeEM.cpp -o BeEM && cd ..

# make the binaries executable
chmod 755 multiprot/multiprot.Linux fiberdock/FiberDock fiberdock/nma
cd ..
```

Check the tool paths in `prism.ini` (`[External_Tools]`) match the unpacked
locations. `naccess` is the default surface tool (`external_tool_choice = 1`).

The full template library (`template/interfaces`, `template/contact`,
`template/hotspot`) and `template_default` ship with this tree — nothing to build.

## Running

```bash
python prism.py <pair_list> <template_list> <jobId>
```

- **`pair_list`** — one target pair per line, `PDBID CHAIN` tokens (≥4 chars each):
  ```
  1cew 2ghuD
  ```
- **`template_list`** — one template interface per line (e.g. `3i6eEF`). Use the
  bundled `template_list` for the default panel.
- **`jobId`** — any name; results land in `jobs/<jobId>/`.

Example:
```bash
python prism.py pair_list template_list mytest
```

## Output

Everything is under `jobs/<jobId>/`:
- `pdb/` — fetched/converted structures
- `surfaceExtract/`, `alignment/`, `transformation/` — stage intermediates
- FiberDock refinement outputs + final interaction energies (also printed to stdout).

Downloaded structures are cached in the top-level `pdb/` folder
(`[Pdb_Folder] pdb_path` in `prism.ini`) and reused across jobs.

## Notes

- The MultiProt (`multiprot.Linux`) and FiberDock binaries are **Linux** builds;
  run on the cluster, not on macOS.
- mmCIF fallback only triggers for entries without a legacy PDB file; it needs the
  compiled `BeEM` binary and network access to `files.pdbj.org`.
