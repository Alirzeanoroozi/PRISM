import sys
import glob

def merge_pdb_files(pdb_files, output_pdb):
    if not pdb_files:
        print("No PDB files provided")
        return

    out = open(output_pdb, 'w')

    first_file = True

    for pdb in pdb_files:
        with open(pdb, 'r') as fh:
            for line in fh:
                record = line[0:6].strip()

                if first_file:
                    # Take everything until MASTER
                    if record == "MASTER":
                        break
                    out.write(line)
                else:
                    # Only ATOM / HETATM lines
                    if record == "ATOM":
                        out.write(line)

        first_file = False

    out.write("END\n")
    out.close()

    print("Merged {} files into {}".format(len(pdb_files), output_pdb))


if __name__ == "__main__":
    # Example usage:
    # python merge_pdb_bundles.py "*_bundle*.pdb" merged.pdb

    if len(sys.argv) != 3:
        print("Usage: python merge_pdb_bundles.py '<glob>' output.pdb")
        sys.exit(1)

    pattern = sys.argv[1]
    output = sys.argv[2]

    pdb_files = sorted(glob.glob(pattern))
    merge_pdb_files(pdb_files, output)
