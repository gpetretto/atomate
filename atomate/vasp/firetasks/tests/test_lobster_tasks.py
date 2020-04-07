import json
import os
import unittest
import shutil
import gridfs
from unittest import mock

from atomate.utils.testing import AtomateTest
from atomate.utils.testing import DB_DIR
from atomate.vasp.firetasks.lobster_tasks import WriteLobsterinputfromIO, LobsterRunToDb, RunLobster
from atomate.vasp.database import VaspCalcDb
from fireworks.utilities.fw_serializers import load_object
from monty.shutil import copy_r
from monty.serialization import dumpfn
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
        cls.vasp_dir = os.path.join(module_dir, "./../../test_files", "lobster", "VASP_calc_for_Lobster")

    def setUp(self):
        super(TestWriteLobsterinputfromIO, self).setUp(lpad=False)

    def test_ioset_explicit(self):
        for fn in ["POSCAR.gz", "POTCAR.gz", "INCAR.gz"]:
            shutil.copy2(os.path.join(self.vasp_dir, fn), ".")
        ft = WriteLobsterinputfromIO(
            poscaraddress="POSCAR.gz", potcaraddress="POTCAR.gz", incaraddress="INCAR.gz", option="standard")
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})

        self.assertEqual(Lobsterin.from_file("lobsterin"), self.ref_lobsterin)

    def test_ioset_settings(self):
        for fn in ["POSCAR.gz", "POTCAR.gz", "INCAR.gz"]:
            shutil.copy2(os.path.join(self.vasp_dir, fn), ".")
        # user supplied lobsterin inputs
        ft = WriteLobsterinputfromIO(
            poscaraddress="POSCAR.gz", potcaraddress="POTCAR.gz", incaraddress="INCAR.gz",
            option="standard",
            user_lobsterin_settings={"COHPEndEnergy": 10.0})
        ft = load_object(ft.to_dict())  # simulate database insertion
        ft.run_task({})

        self.assertEqual(Lobsterin.from_file("lobsterin"), self.ref_lobsterin2)


class TestLobsterRunToDb(AtomateTest):
    @classmethod
    def setUpClass(cls):
        cls.vasp_dir = os.path.join(module_dir, "./../../test_files", "lobster", "vasp_lobster_output")
        cls.vasp_si_dir = os.path.join(module_dir, "./../../test_files", "lobster", "si_vasp_lobster2/lobster/outputs")

    def setUp(self):
        super(TestLobsterRunToDb, self).setUp(lpad=False)

    def test_jsonfile(self):
        copy_r(self.vasp_dir, ".")
        ft = LobsterRunToDb(calc_loc=True)
        ft.run_task({"calc_locs": [{"path": "test"}]})
        with open("task_lobster.json") as f:
            load_dict = json.load(f)
        self.assertEqual(load_dict["formula_pretty"], 'K2Sn2O3')
        self.assertListEqual(load_dict["output"]["chargespilling"], [0.008, 0.008])

    def test_mongodb(self):
        try:
            VaspCalcDb.from_db_file(DB_FILE)
        except:
            raise unittest.SkipTest('Cannot connect to MongoDB! Is the database server running? '
                                    'Are the credentials correct?')

        copy_r(self.vasp_dir, ".")
        ft = LobsterRunToDb(calc_loc=True, db_file=DB_FILE, store_all_outputs=False)
        ft.run_task({"calc_locs": [{"path": "test"}]})
        coll = self.get_task_collection("lobster")
        load_dict = coll.find_one({"formula_pretty": "K2Sn2O3"})
        self.assertEqual(load_dict["formula_pretty"], 'K2Sn2O3')
        self.assertListEqual(load_dict["output"]["chargespilling"], [0.008, 0.008])
        self.assertNotIn("lobster_icohplist", load_dict)

    def test_mongodb_all_files(self):
        try:
            VaspCalcDb.from_db_file(DB_FILE)
        except:
            raise unittest.SkipTest('Cannot connect to MongoDB! Is the database server running? '
                                    'Are the credentials correct?')

        copy_r(self.vasp_dir, ".")
        ft = LobsterRunToDb(calc_loc=True, db_file=DB_FILE, store_all_outputs=True)
        ft.run_task({"calc_locs": [{"path": "test"}]})
        coll = self.get_task_collection("lobster")
        load_dict = coll.find_one({"formula_pretty": "K2Sn2O3"})
        self.assertEqual(load_dict["formula_pretty"], 'K2Sn2O3')
        self.assertListEqual(load_dict["output"]["chargespilling"], [0.008, 0.008])
        db = self.get_task_database()
        for fn in ["ICOHPLIST", "ICOOPLIST", "COHPCAR",
                   "COOPCAR", "GROSSPOP", "CHARGE", "DOSCAR"]:
            oid = load_dict[fn.lower() + "_id"]
            gfs = gridfs.GridFS(db, "lobster_" + fn.lower())
            results = gfs.find({"_id": oid}).count()
            self.assertEqual(results, 1)

    def test_jsonfile_si(self):
        copy_r(self.vasp_si_dir, ".")
        ft = LobsterRunToDb(calc_loc=True)
        ft.run_task({"calc_locs": [{"path": "test"}]})
        with open("task_lobster.json") as f:
            load_dict = json.load(f)
        self.assertEqual(load_dict["formula_pretty"], 'Si')
        self.assertListEqual(load_dict["output"]["chargespilling"], [0.0141, 0.0141])


class TestRunLobster(AtomateTest):

    def setUp(self):
        super(TestRunLobster, self).setUp(lpad=False)

    def test_run_task(self):
        with mock.patch("atomate.vasp.firetasks.lobster_tasks.Custodian") as mock_c:
            instance = mock_c.return_value
            instance.run.return_value = 'test'
            t = RunLobster(lobster_cmd='', gzip_output=True, gzip_WAVECAR=False)
            self.assertIsNone(t.run_task(fw_spec={}))
            mock_c.assert_called_once()
            self.assertEqual(len(mock_c.call_args[1]["handlers"]), 0)
            self.assertEqual(len(mock_c.call_args[1]["validators"]), 2)

            mock_c.reset_mock()
            t = RunLobster(lobster_cmd='', gzip_output=False, gzip_WAVECAR=False, handler_group="no_handler",
                           validator_group="strict")
            self.assertIsNone(t.run_task(fw_spec={}))
            mock_c.assert_called_once()
            self.assertEqual(len(mock_c.call_args[1]["handlers"]), 0)
            self.assertEqual(len(mock_c.call_args[1]["validators"]), 3)

            mock_c.reset_mock()
            t = RunLobster(lobster_cmd='', gzip_output=False, gzip_WAVECAR=True, handler_group=[],
                           validator_group="no_validator")
            self.assertIsNone(t.run_task(fw_spec={}))
            mock_c.assert_called_once()
            self.assertEqual(len(mock_c.call_args[1]["handlers"]), 0)
            self.assertEqual(len(mock_c.call_args[1]["validators"]), 0)

            mock_c.reset_mock()
            d = {"test": True}
            dumpfn(d, "custodian.json")
            t = RunLobster(lobster_cmd='', gzip_output=False, gzip_WAVECAR=True)
            a = t.run_task(fw_spec={})
            self.assertDictEqual(a.stored_data["custodian"], d)
            mock_c.assert_called_once()


if __name__ == '__main__':
    unittest.main()
