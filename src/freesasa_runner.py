"""Run FreeSASA in a configured Python environment and emit residue areas."""

import json
import sys

import freesasa


def main(argv):
    if len(argv) != 3:
        raise SystemExit("usage: python -m src.freesasa_runner INPUT.pdb OUTPUT.json")

    input_path, output_path = argv[1:]
    structure = freesasa.Structure(input_path)
    parameters = freesasa.Parameters(
        {"algorithm": "LeeRichards", "probe-radius": 1.4, "n-slices": 20}
    )
    result = freesasa.calc(structure, parameters)
    areas = {}
    for chain, residues in result.residueAreas().items():
        areas[chain] = {}
        for residue_number, residue in residues.items():
            areas[chain][residue_number] = {
                "residue_name": residue.residueType,
                "total": residue.total,
            }

    with open(output_path, "w") as handle:
        json.dump(areas, handle, sort_keys=True)


if __name__ == "__main__":
    main(sys.argv)
