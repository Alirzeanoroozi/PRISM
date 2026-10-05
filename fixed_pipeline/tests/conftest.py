"""Configure PYTHONPATH for isolated module imports."""
import sys
from pathlib import Path

# Add ORIGINAL src/ directory FIRST so that shared dependencies
# (contact, candidate_audit, pdb_download, utils, etc.) resolve.
# The fixed modules in fixed_pipeline/src/ use relative imports like
# `from .contact import ...` — those resolve via the `src` package.
# But contact, candidate_audit, etc. are only in the ORIGINAL src/.
orig_src = str(Path(__file__).resolve().parent.parent.parent / "src")
if orig_src not in sys.path:
    sys.path.insert(0, orig_src)

# Add fixed_pipeline/src so fixed modules override originals.
fixed_src = str(Path(__file__).resolve().parent.parent / "src")
if fixed_src not in sys.path:
    sys.path.insert(0, fixed_src)
