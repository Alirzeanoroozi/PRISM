import argparse
import csv
import os


def parse_ply_vertices(path):
    with open(path, "r") as handle:
        line = handle.readline().strip()
        if line != "ply":
            raise ValueError("Only ASCII PLY files are supported")

        vertex_count = None
        properties = []
        in_vertex = False

        for line in handle:
            stripped = line.strip()
            if stripped.startswith("element "):
                parts = stripped.split()
                in_vertex = len(parts) >= 3 and parts[1] == "vertex"
                if in_vertex:
                    vertex_count = int(parts[2])
                    properties = []
            elif stripped.startswith("property ") and in_vertex:
                parts = stripped.split()
                properties.append(parts[-1])
            elif stripped == "end_header":
                break

        if vertex_count is None:
            raise ValueError("No vertex element found in PLY header")
        if not {"x", "y", "z"}.issubset(set(properties)):
            raise ValueError("PLY vertex properties must include x, y, z")

        rows = []
        for _ in range(vertex_count):
            values = handle.readline().strip().split()
            if len(values) != len(properties):
                raise ValueError("PLY vertex row does not match header property count")
            row = dict(zip(properties, values))
            rows.append({"x": row["x"], "y": row["y"], "z": row["z"]})
        return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ply", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = parse_ply_vertices(args.ply)
    parent = os.path.dirname(args.output)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["x", "y", "z"])
        writer.writeheader()
        writer.writerows(rows)
    print(args.output)


if __name__ == "__main__":
    main()
