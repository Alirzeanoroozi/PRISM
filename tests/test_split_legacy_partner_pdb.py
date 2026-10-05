from benchmark.scripts.split_legacy_partner_pdb import split_pdb


def _atom(serial, residue):
    return (
        f"ATOM  {serial:5d}  CA  ALA B{residue:4d}    "
        "0.000   0.000   0.000  1.00 20.00           C\n"
    )


def test_split_legacy_partner_pdb_records_boundary_and_distinct_chains(tmp_path):
    source = tmp_path / "legacy.pdb"
    destination = tmp_path / "repaired.pdb"
    source.write_text(_atom(1, 1) + _atom(2, 2) + _atom(3, 1) + "END\n")
    record = split_pdb(source, destination)
    assert record["first_partner_atom_records"] == 2
    assert record["second_partner_atom_records"] == 1
    text = destination.read_text()
    assert text.splitlines()[0][21] == "A"
    assert text.splitlines()[2][21] == "B"
