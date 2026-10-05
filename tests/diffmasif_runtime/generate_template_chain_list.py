import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--templates", default="templates/templates.txt")
    parser.add_argument("--output", default="tests/diffmasif_runtime/template_chain_list.txt")
    args = parser.parse_args()

    template_ids = [
        line.strip()
        for line in Path(args.templates).read_text().splitlines()
        if line.strip()
    ]
    chains = set()
    for template_id in template_ids:
        pdb_id = template_id[:4].upper()
        chains.add(f"{pdb_id}_{template_id[4]}")
        chains.add(f"{pdb_id}_{template_id[5]}")

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(sorted(chains)) + "\n")
    print(output_path)
    print(len(chains))


if __name__ == "__main__":
    main()
