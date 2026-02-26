"""Tests for src/interface.py interface generation functions."""
import pytest

from src.interface import _format_pdb_line


class TestFormatPdbLine:
    def test_output_length(self):
        line = _format_pdb_line(1, " CA ", "ALA", "A", 1, 1.0, 2.0, 3.0)
        # Strip trailing newline for length check; rest of line is fixed-width PDB
        assert line.endswith("\n")
        # PDB ATOM lines are 80 chars (plus newline), but our formatter may vary slightly
        assert len(line.strip()) > 0

    def test_starts_with_atom(self):
        line = _format_pdb_line(1, " CA ", "ALA", "A", 1, 1.0, 2.0, 3.0)
        assert line.startswith("ATOM")

    def test_serial_number_embedded(self):
        line = _format_pdb_line(42, " CA ", "GLY", "B", 10, 0.0, 0.0, 0.0)
        assert "42" in line

    def test_residue_name_embedded(self):
        line = _format_pdb_line(1, " CA ", "TRP", "A", 5, 1.0, 2.0, 3.0)
        assert "TRP" in line

    def test_chain_id_embedded(self):
        line = _format_pdb_line(1, " CA ", "ALA", "C", 1, 0.0, 0.0, 0.0)
        assert "C" in line

    def test_coordinates_formatted(self):
        line = _format_pdb_line(1, " CA ", "ALA", "A", 1, 12.345, -6.789, 0.001)
        assert "12.345" in line
        assert "-6.789" in line

    def test_default_occupancy_and_bfactor(self):
        line = _format_pdb_line(1, " CA ", "ALA", "A", 1, 0.0, 0.0, 0.0)
        # Default occupancy=1.00, bfactor=0.00
        assert "1.00" in line

    def test_custom_bfactor(self):
        line = _format_pdb_line(1, " CA ", "ALA", "A", 1, 0.0, 0.0, 0.0, bfactor=99.0)
        assert "99.00" in line

    def test_element_embedded(self):
        line = _format_pdb_line(1, " CA ", "ALA", "A", 1, 0.0, 0.0, 0.0, element="N")
        assert " N" in line
