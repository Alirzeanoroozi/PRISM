"""Tests for TM-score parsing in benchmark/scripts/test_iptm.py."""
import sys
import os
import pytest

# Add benchmark/scripts to path so we can import the module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "benchmark", "scripts"))
from test_iptm import _parse_tmalign_output


# Sample TMalign output (condensed real output format).
# Alignment statistics ("Aligned length", "TM-score") appear at column 0
# in actual TMalign output, so the sample reflects that exact format.
_SAMPLE_OUTPUT = """\
 ******************************************************************************
 * TM-align (Version 20190822): A protein structural alignment algorithm      *
 ******************************************************************************

 Name of Chain_1: input.pdb (to be superimposed onto Chain_2)
 Name of Chain_2: reference.pdb

Aligned length= 51, RMSD=   0.00, Seq_ID=n_identical/n_aligned= 1.000
TM-score= 0.85321 (if normalized by length of Chain_1, i.e., L=51, d0=3.25)
TM-score= 0.90132 (if normalized by length of Chain_2, i.e., L=55, d0=3.38)

 (":"  denotes residue pairs of distance < 5.0 Angstrom)
ACDEF
:::::
ACDEF
"""

_SINGLE_CHAIN_OUTPUT = """\
Aligned length= 10, RMSD=   1.23, Seq_ID=n_identical/n_aligned= 0.750
TM-score= 0.72000 (if normalized by length of Chain_1, i.e., L=10, d0=2.0)
TM-score= 0.68000 (if normalized by length of Chain_2, i.e., L=12, d0=2.1)
"""


class TestParseTmalignOutput:
    def test_parses_aligned_length(self):
        result = _parse_tmalign_output(_SAMPLE_OUTPUT)
        assert result["aligned_length"] == 51

    def test_parses_rmsd(self):
        result = _parse_tmalign_output(_SAMPLE_OUTPUT)
        assert result["rmsd"] == pytest.approx(0.00)

    def test_parses_seq_id(self):
        result = _parse_tmalign_output(_SAMPLE_OUTPUT)
        assert result["seq_id"] == pytest.approx(1.000)

    def test_parses_tmscore_query(self):
        result = _parse_tmalign_output(_SAMPLE_OUTPUT)
        assert result["tmscore_query"] == pytest.approx(0.85321)

    def test_parses_tmscore_ref(self):
        result = _parse_tmalign_output(_SAMPLE_OUTPUT)
        assert result["tmscore_ref"] == pytest.approx(0.90132)

    def test_empty_output_returns_nones(self):
        result = _parse_tmalign_output("")
        assert result["tmscore_ref"] is None
        assert result["tmscore_query"] is None
        assert result["rmsd"] is None
        assert result["aligned_length"] is None

    def test_single_chain_parses_scores(self):
        result = _parse_tmalign_output(_SINGLE_CHAIN_OUTPUT)
        assert result["tmscore_query"] == pytest.approx(0.72)
        assert result["tmscore_ref"] == pytest.approx(0.68)
        assert result["rmsd"] == pytest.approx(1.23)

    def test_return_keys_present(self):
        result = _parse_tmalign_output(_SAMPLE_OUTPUT)
        for key in ("tmscore_ref", "tmscore_query", "rmsd", "aligned_length", "seq_id"):
            assert key in result
