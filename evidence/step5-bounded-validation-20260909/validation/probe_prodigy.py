import json
import pathlib
import sys

from src import prodigy_ranker


left, right, output_dir = sys.argv[1:]
selected = prodigy_ranker.select_top_candidates_with_prodigy(
    [(left, right)],
    top_k=1,
    executable="prodigy",
    output_dir=output_dir,
)
state_path = pathlib.Path(output_dir) / "prodigy-state.jsonl"
states = [
    json.loads(line)["state"]
    for line in state_path.read_text(encoding="utf-8").splitlines()
    if line.strip()
]
score_files = sorted(pathlib.Path(output_dir).glob("*.json"))
score_payload = json.loads(score_files[0].read_text()) if score_files else {}
print(json.dumps({
    "selected_count": len(selected),
    "states": states,
    "score_status": score_payload.get("status"),
    "affinity_kcal_mol": score_payload.get("affinity_kcal_mol"),
}, sort_keys=True))
