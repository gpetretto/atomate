import json
import os
import unittest
from unittest import mock

from atomate.utils.testing import AtomateTest
from atomate.utils.testing import DB_DIR
from atomate.vasp.firetasks.lobster_tasks import WriteLobsterinputfromIO, LobsterRunToDb, RunLobster
from fireworks.utilities.fw_serializers import load_object
from monty.os import cd
from monty.tempfile import ScratchDir
from pymatgen.io.lobster import Lobsterin

DB_FILE = os.path.join(DB_DIR, "db.json")
module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class TestWriteLobsterinputfromIO(AtomateTest):
    @classmethod
    def setUpClass(cls):
        cls.ref_lobsterin = Lobsterin.from_file(
            os.path.join(module_dir, "./../../test_files", "lobster", "Lobsterinputs", "lobsterin"))
        cls.ref_lobsterin2 = Lobsterin.from_file(
            os.path.join(module_dir, "./../../test_files", "lobster", "lobsterins", "lobsterin2"))

    def setUp(self):
        super(TestWriteLobsterinputfromIO, self).setUp(lpad=False)

    def test_ioset_explicit(self):
        with cd(os.path.join(module_dir, "./../../test_files", "lobster", "VASP_calc_for_Lobster")):
            with ScratchDir('.', copy_from_current_on_enter=True) as d:
                ft = WriteLobsterinputfromIO(
                    poscaraddress="POSCAR.gz", potcaraddress="POTCAR.gz", incaraddress="INCAR.gz", option="standard")
                ft = load_object(ft.to_dict())  # simulate database insertion
                ft.run_task({})

                self.assertEqual(Lobsterin.from_file("lobsterin"), self.ref_lobsterin)
            with ScratchDir('.', copy_from_current_on_enter=True) as d:
                # user supplied lobsterin inputs
                ft = WriteLobsterinputfromIO(
                    poscaraddress="POSCAR.gz", potcaraddress="POTCAR.gz", incaraddress="INCAR.gz",
                    option="standard",
                    user_lobsterin_settings={"COHPEndEnergy": 10.0})
                ft = load_object(ft.to_dict())  # simulate database insertion
                ft.run_task({})

                self.assertEqual(Lobsterin.from_file("lobsterin"), self.ref_lobsterin2)

    def tearDown(self):
        pass


class TestLobsterRunToDb(AtomateTest):

    def setUp(self):
        pass

    def test_generatedbfile(self):
        with cd(os.path.join(module_dir, "./../../test_files", "lobster", "vasp_lobster_output")):
            with ScratchDir('.', copy_from_current_on_enter=True) as d:
                ft = LobsterRunToDb(calc_loc=True)
                ft.run_task({"calc_locs": [{"path": "test"}]})
                with open("task_lobster.json") as f:
                    load_dict = json.load(f)
                self.assertEqual(load_dict["formula_pretty"], 'K2Sn2O3')
                self.assertListEqual(load_dict["output"]["chargespilling"], [0.008, 0.008])

        with cd(os.path.join(module_dir, "./../../test_files", "lobster", "vasp_lobster_output")):
            with ScratchDir('.', copy_from_current_on_enter=True) as d:
                ft = LobsterRunToDb(calc_loc=True, db_file=DB_FILE)
                ft.run_task({"calc_locs": [{"path": "test"}]})
                coll = self.get_task_collection("lobster")
                load_dict = coll.find_one({"formula_pretty": "K2Sn2O3"})
                # will this always return the correct entry?
                self.assertEqual(load_dict["formula_pretty"], 'K2Sn2O3')
                self.assertListEqual(load_dict["output"]["chargespilling"], [0.008, 0.008])

        with cd(os.path.join(module_dir, "./../../test_files", "lobster", "si_vasp_lobster2/lobster/outputs")):
            with ScratchDir('.', copy_from_current_on_enter=True) as d:
                ft = LobsterRunToDb(calc_loc=True)
                ft.run_task({"calc_locs": [{"path": "test"}]})
                with open("task_lobster.json") as f:
                    load_dict = json.load(f)
                self.assertEqual(load_dict["formula_pretty"], 'Si')
                self.assertListEqual(load_dict["output"]["chargespilling"], [0.0141, 0.0141])

    def tearDown(self):
        db = self.get_task_database()
        for coll in db.collection_names():
            if coll != "system.indexes":
                db[coll].drop()


class TestRunLobster(unittest.TestCase):
    def setUp(self) -> None:
        pass

    def test_run_task(self):
        with mock.patch("atomate.vasp.firetasks.lobster_tasks.Custodian") as mock_c:
            instance = mock_c.return_value
            instance.run_task.return_value = 'test'
            t = RunLobster(lobster_cmd='', gzip_output=True, gzip_WAVECAR=False, strict_handlers_validators=False)
            t.run_task(fw_spec={})
            t = RunLobster(lobster_cmd='', gzip_output=False, gzip_WAVECAR=False, strict_handlers_validators=False)
            t.run_task(fw_spec={})
            t = RunLobster(lobster_cmd='', gzip_output=False, gzip_WAVECAR=True, strict_handlers_validators=False)
            t.run_task(fw_spec={})
            t = RunLobster(lobster_cmd='', gzip_output=False, gzip_WAVECAR=True, strict_handlers_validators=True)
            t.run_task(fw_spec={})



if __name__ == '__main__':
    unittest.main()
