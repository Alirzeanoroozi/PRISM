import argparse
import csv
import json
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

STANDARD_ASA = {
    "ALA": 107.95,
    "CYS": 134.28,
    "ASP": 140.39,
    "GLU": 172.25,
    "PHE": 199.48,
    "GLY": 80.10,
    "HIS": 182.88,
    "ILE": 175.12,
    "LYS": 200.81,
    "LEU": 178.63,
    "MET": 194.15,
    "ASN": 143.94,
    "PRO": 136.13,
    "GLN": 178.50,
    "ARG": 238.76,
    "SER": 116.50,
    "THR": 139.27,
    "VAL": 151.44,
    "TRP": 249.36,
    "TYR": 212.76,
}


def _infer_element(atom_name: str, element_field: str) -> str:
    element = (element_field or "").strip()
    if element:
        return element.upper()
    atom_name = atom_name.strip()
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    return (letters[:1] or "").upper()


def parse_pdb_residues(pdb_path: str, allowed_chains: Optional[Iterable[str]] = None) -> List[Dict[str, object]]:
    chain_filter = set(allowed_chains or [])
    grouped: Dict[Tuple[str, int, str], List[Tuple[float, float, float]]] = {}

    with open(pdb_path, "r") as handle:
        for line in handle:
            record = line[:6].strip()
            if record != "ATOM":
                continue
            chain_id = line[21].strip()
            if chain_filter and chain_id not in chain_filter:
                continue
            resname = line[17:20].strip().upper()
            if resname not in STANDARD_AA:
                continue
            altloc = line[16].strip()
            if altloc not in ("", "A"):
                continue
            atom_name = line[12:16].strip()
            element = _infer_element(atom_name, line[76:78])
            if element == "H":
                continue
            resseq = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            grouped.setdefault((chain_id, resseq, resname), []).append((x, y, z))

    residues: List[Dict[str, object]] = []
    for (chain_id, resseq, resname), atoms in grouped.items():
        residues.append(
            {
                "key": f"{resname}_{resseq}_{chain_id}",
                "resname": resname,
                "resseq": resseq,
                "chain_id": chain_id,
                "atom_coords": atoms,
            }
        )
    residues.sort(key=lambda item: (item["chain_id"], item["resseq"], item["resname"]))
    return residues


def load_surface_points_csv(path: str, score_column: str = "score") -> List[Dict[str, float]]:
    points: List[Dict[str, float]] = []
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = {"x", "y", "z"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Surface point CSV is missing required columns: {sorted(missing)}")
        for row in reader:
            point = {
                "x": float(row["x"]),
                "y": float(row["y"]),
                "z": float(row["z"]),
            }
            if score_column in row and row[score_column] not in ("", None):
                point["score"] = float(row[score_column])
            points.append(point)
    if not points:
        raise ValueError(f"No surface points were loaded from {path}")
    return points


def squared_distance(a: Sequence[float], b: Sequence[float]) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def nearest_residue(
    point: Sequence[float],
    residues: Sequence[Dict[str, object]],
    assignment_cutoff: float,
) -> Optional[Dict[str, object]]:
    best_residue: Optional[Dict[str, object]] = None
    best_distance_sq = assignment_cutoff ** 2
    for residue in residues:
        for atom_coord in residue["atom_coords"]:
            dist_sq = squared_distance(point, atom_coord)
            if dist_sq <= best_distance_sq:
                best_distance_sq = dist_sq
                best_residue = residue
    return best_residue


def build_proxy_accessibility(
    pdb_path: str,
    surface_points_path: str,
    allowed_chains: Optional[Iterable[str]] = None,
    assignment_cutoff: float = 4.0,
    score_column: str = "score",
) -> Dict[str, Dict[str, float]]:
    residues = parse_pdb_residues(pdb_path, allowed_chains=allowed_chains)
    if not residues:
        raise ValueError(f"No standard amino-acid residues found in {pdb_path}")

    points = load_surface_points_csv(surface_points_path, score_column=score_column)
    return build_proxy_accessibility_from_points(
        pdb_path=pdb_path,
        points=points,
        allowed_chains=allowed_chains,
        assignment_cutoff=assignment_cutoff,
    )


def build_proxy_accessibility_from_points(
    pdb_path: str,
    points: Sequence[Dict[str, float]],
    allowed_chains: Optional[Iterable[str]] = None,
    assignment_cutoff: float = 4.0,
) -> Dict[str, Dict[str, float]]:
    residues = parse_pdb_residues(pdb_path, allowed_chains=allowed_chains)
    if not residues:
        raise ValueError(f"No standard amino-acid residues found in {pdb_path}")
    if not points:
        raise ValueError("No surface points were provided")

    residue_stats: Dict[str, Dict[str, float]] = {
        residue["key"]: {
            "point_count": 0.0,
            "score_sum": 0.0,
            "score_count": 0.0,
            "point_density_per_standard_asa": 0.0,
            "relative_accessibility_proxy": 0.0,
        }
        for residue in residues
    }

    unassigned_points = 0
    for point in points:
        residue = nearest_residue(
            (point["x"], point["y"], point["z"]),
            residues,
            assignment_cutoff=assignment_cutoff,
        )
        if residue is None:
            unassigned_points += 1
            continue
        stats = residue_stats[residue["key"]]
        stats["point_count"] += 1.0
        if "score" in point:
            stats["score_sum"] += point["score"]
            stats["score_count"] += 1.0

    max_density = 0.0
    for residue in residues:
        stats = residue_stats[residue["key"]]
        standard_asa = STANDARD_ASA.get(residue["resname"])
        if standard_asa:
            stats["point_density_per_standard_asa"] = stats["point_count"] / standard_asa
            max_density = max(max_density, stats["point_density_per_standard_asa"])

    if max_density > 0.0:
        for stats in residue_stats.values():
            stats["relative_accessibility_proxy"] = (
                100.0 * stats["point_density_per_standard_asa"] / max_density
            )

    for stats in residue_stats.values():
        stats["mean_score"] = (
            stats["score_sum"] / stats["score_count"] if stats["score_count"] > 0 else None
        )

    residue_stats["_meta"] = {
        "assignment_cutoff": assignment_cutoff,
        "surface_point_count": float(len(points)),
        "assigned_surface_point_count": float(len(points) - unassigned_points),
        "unassigned_surface_point_count": float(unassigned_points),
    }
    return residue_stats


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True, help="Path to the structure PDB file.")
    parser.add_argument("--surface-points", required=True, help="CSV file with x,y,z point columns.")
    parser.add_argument("--output", required=True, help="Output JSON path.")
    parser.add_argument("--chains", nargs="*", default=None)
    parser.add_argument("--assignment-cutoff", type=float, default=4.0)
    parser.add_argument("--score-column", default="score")
    args = parser.parse_args()

    result = build_proxy_accessibility(
        pdb_path=args.pdb,
        surface_points_path=args.surface_points,
        allowed_chains=args.chains,
        assignment_cutoff=args.assignment_cutoff,
        score_column=args.score_column,
    )
    with open(args.output, "w") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
