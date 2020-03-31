# coding: utf-8


import os
import unittest

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import RotationTransformation

from atomate.vasp.fireworks.core import StaticFW
from atomate.vasp.fireworks.lobster import LobsterFW

__author__ = 'Janine George, Guido Petretto'
__email__ = 'janine.george@uclouvain.be'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_dir = os.path.join(module_dir, "..", "..", "..", "common", "test_files")
reference_dir = os.path.join(module_dir, "..", "..", "test_files")

DEBUG_MODE = False  # If true, retains the database and output dirs at the end of the test
VASP_CMD = None  # If None, runs a "fake" VASP. Otherwise, runs VASP with this command...


class TestLobsterFireworks(unittest.TestCase):
    def setUp(self):
        coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
        lattice = Lattice([[3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00],
                           [0.00, -2.2171384943, 3.1355090603]])
        self.structure = Structure(lattice, ["Si", "Si"], coords)

    def testStaticFW(self):
        static_fw=StaticFW(structure=self.structure).name
        self.assertEqual(LobsterFW(structure=self.structure, parents=static_fw).name, "Si-lobster_calculation")
        lobster_fw=LobsterFW(prev_calc_dir="/")
        self.assertEqual(lobster_fw.name, "unknown-lobster_calculation")
        self.assertEqual(lobster_fw.tasks[0]["calc_dir"],"/")
        #check for ValueError when no parent or calc_dir are provided
        with self.assertRaises(ValueError):
            LobsterFW()
 













if __name__ == "__main__":
    unittest.main()
