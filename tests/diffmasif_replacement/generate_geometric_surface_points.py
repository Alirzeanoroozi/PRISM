import argparse
import csv
import math
from typing import Iterable, List, Optional, Tuple


STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

VDW_RADIUS = {
    "C": 1.76,
    "N": 1.65,
    "O": 1.40,
    "S": 1.85,
    "P": 1.90,
}

SAMPLE_DIRECTIONS = [
    (1.0, 0.0, 0.0),
    (-1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, -1.0, 0.0),
    (0.0, 0.0, 1.0),
    (0.0, 0.0, -1.0),
    (1.0, 1.0, 1.0),
    (1.0, 1.0, -1.0),
    (1.0, -1.0, 1.0),
    (1.0, -1.0, -1.0),
    (-1.0, 1.0, 1.0),
    (-1.0, 1.0, -1.0),
    (-1.0, -1.0, 1.0),
    (-1.0, -1.0, -1.0),
]


def _infer_element(atom_name: str, element_field: str) -> str:
    element = (element_field or "").strip()
    if element:
        return element.upper()
    letters = "".join(ch for ch in atom_name.strip() if ch.isalpha())
    return (letters[:1] or "").upper()


def parse_atoms(pdb_path: str, allowed_chains: Optional[Iterable[str]] = None) -> List[dict]:
    allowed = set(allowed_chains or [])
    atoms: List[dict] = []
    with open(pdb_path, "r") as handle:
        for line in handle:
            if line[:6].strip() != "ATOM":
                continue
            chain_id = line[21].strip()
            if allowed and chain_id not in allowed:
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
            atoms.append(
                {
                    "chain_id": chain_id,
                    "resname": resname,
                    "resseq": int(line[22:26].strip()),
                    "element": element,
                    "x": float(line[30:38].strip()),
                    "y": float(line[38:46].strip()),
                    "z": float(line[46:54].strip()),
                }
            )
    return atoms


def _normalize(direction: Tuple[float, float, float]) -> Tuple[float, float, float]:
    norm = math.sqrt(direction[0] ** 2 + direction[1] ** 2 + direction[2] ** 2)
    return (direction[0] / norm, direction[1] / norm, direction[2] / norm)


NORMALIZED_DIRECTIONS = [_normalize(direction) for direction in SAMPLE_DIRECTIONS]


def _dist_sq(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def generate_surface_points(
    pdb_path: str,
    output_csv: str,
    allowed_chains: Optional[Iterable[str]] = None,
    probe_radius: float = 1.4,
    occlusion_tolerance: float = 0.2,
) -> int:
    atoms = parse_atoms(pdb_path, allowed_chains=allowed_chains)
    rows: List[Tuple[float, float, float, str, int, str]] = []

    atom_data = []
    for atom in atoms:
        radius = VDW_RADIUS.get(atom["element"], 1.8)
        atom_data.append((atom, radius))

    for atom, radius in atom_data:
        shell_radius = radius + probe_radius
        for direction in NORMALIZED_DIRECTIONS:
            point = (
                atom["x"] + shell_radius * direction[0],
                atom["y"] + shell_radius * direction[1],
                atom["z"] + shell_radius * direction[2],
            )
            occluded = False
            for other_atom, other_radius in atom_data:
                if other_atom is atom:
                    continue
                limit = max(other_radius + probe_radius - occlusion_tolerance, 0.1)
                if _dist_sq(point, (other_atom["x"], other_atom["y"], other_atom["z"])) < limit ** 2:
                    occluded = True
                    break
            if not occluded:
                rows.append((point[0], point[1], point[2], atom["chain_id"], atom["resseq"], atom["resname"]))

    with open(output_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x", "y", "z", "chain", "resseq", "resname"])
        writer.writerows(rows)
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--chains", nargs="*", default=None)
    parser.add_argument("--probe-radius", type=float, default=1.4)
    parser.add_argument("--occlusion-tolerance", type=float, default=0.2)
    args = parser.parse_args()

    count = generate_surface_points(
        pdb_path=args.pdb,
        output_csv=args.output,
        allowed_chains=args.chains,
        probe_radius=args.probe_radius,
        occlusion_tolerance=args.occlusion_tolerance,
    )
    print(count)


if __name__ == "__main__":
    main()
