#!/usr/bin/env python2.7
import os
import shutil
import tempfile
import unittest
import new

from structuralAlignmentTM import StructuralAligner


class StructuralAlignmentTMTests(unittest.TestCase):
    def test_parse_tmalign_ignores_blank_matrix_lines(self):
        original_cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp()
        try:
            os.chdir(tmpdir)
            os.mkdir("alignment")
            with open("matrix.out", "w") as matrix_file:
                matrix_file.write(
                    "------ matrix ------\n"
                    "m t u0 u1 u2\n"
                    "0 1.0 1.0 0.0 0.0\n"
                    "1 2.0 0.0 1.0 0.0\n"
                    "2 3.0 0.0 0.0 1.0\n"
                    "\n"
                    "Code for rotating Structure A\n"
                )
            with open("alignment/out.tm", "w") as tm_file:
                tm_file.write(
                    "Aligned length= 1, RMSD= 0.00, Seq_ID=n_identical/n_aligned= 1.000\n"
                    "TM-score= 0.50000 (if normalized by length of Chain_1, i.e., LN=1, d0=0.50)\n"
                    "TM-score= 0.25000 (if normalized by length of Chain_2, i.e., LN=2, d0=0.50)\n"
                )

            aligner = new.instance(StructuralAligner)
            parsed = aligner.parseTMalign("query.pdb", "interface.pdb")
            self.assertEqual(parsed[0][0], 1)
            self.assertEqual(parsed[0][2][0], [1.0, 2.0, 3.0])
            self.assertEqual(
                parsed[0][2][1],
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            )
            self.assertEqual(parsed[0][4], 0.5)
        finally:
            os.chdir(original_cwd)
            shutil.rmtree(tmpdir)


if __name__ == "__main__":
    unittest.main()
