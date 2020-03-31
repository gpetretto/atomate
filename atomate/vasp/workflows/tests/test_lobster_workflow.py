__author__ = "Janine George, Guido Petretto"
__email__ = "janine.george@uclouvain.be"

import json
import os
import unittest

from atomate.utils.testing import AtomateTest
from atomate.utils.testing import DB_DIR
from atomate.vasp.powerups import use_fake_vasp, use_custodian
from atomate.vasp.workflows.base.lobster import get_all_possible_basis_combinations, get_wf_lobster
from atomate.vasp.firetasks.lobster_tasks import use_fake_lobster
from fireworks.core.rocket_launcher import rapidfire
from pymatgen.util.testing import PymatgenTest

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

VASP_CMD = None
LOBSTER_CMD = None
DB_FILE = os.path.join(DB_DIR, "db.json")
refs_dirs_si_vasp = {"static": os.path.join(module_dir, "../../test_files/lobster/si_vasp_lobster/vasp/")}
refs_dirs_si_vasp2 = {"static": os.path.join(module_dir, "../../test_files/lobster/si_vasp_lobster2/vasp/")}
refs_dirs_si_lobster = {"lobster_calculation": os.path.join(module_dir, "../../test_files/lobster/si_vasp_lobster/lobster")}
refs_dirs_si_lobster2 = {
    "lobster_calculation": os.path.join(module_dir, "../../lobster/test_files/lobster/si_vasp_lobster2/lobster")}



#TODO: work on this part! check that fake_lobster works correctly
class TestWF_Lobster(AtomateTest):
    def setUp(self):
        super(TestWF_Lobster, self).setUp(lpad=True)
        self.struct_si = PymatgenTest.get_structure("Si")

    def _check_run(self, d, mode):
        if mode not in ["lobsternormal"]:
            raise ValueError("Invalid mode!")

        self.assertEqual(d["formula_pretty"], "Si")
        self.assertEqual(d["formula_anonymous"], "A")
        self.assertEqual(d["nelements"], 1)
        self.assertEqual(d["state"], "successful")

        if mode in ["lobsternormal"]:
            self.assertListEqual(d["output"]["chargespilling"], [0.0147, 0.0147])
            self.assertListEqual(d["output"]["elements"], ['Si'])
            self.assertListEqual(d["output"]["basistype"], ['pbeVaspFit2015'])
            self.assertListEqual(d["output"]["basisfunctions"], [['3s', '3p_y', '3p_z', '3p_x']])
            self.assertTrue(d["output"]["hasDOSCAR"])

    def _single_Vasp_Lobster(self, delete_wavecars=False, gzipped_output=False, gzipped_WAVECAR=False, fw_id=1,
                             user_supplied_basis=None, backup=False, fake=True):
        # add the workflow

        structure = self.struct_si
        # with ScratchDir('.', copy_from_current_on_enter=True) as d:
        my_wf = get_wf_lobster(structure=structure, c={"vasp_cmd": VASP_CMD, "DB_FILE": None},
                               user_kpoints_settings={"grid_density": 100}, delete_all_wavecars=delete_wavecars,
                               user_supplied_basis=user_supplied_basis)
        if fake:
            my_wf = use_fake_vasp(my_wf, refs_dirs_si_vasp)
            my_wf = use_fake_lobster(my_wf, refs_dirs_si_lobster, gzipped_output=gzipped_output,
                                     gzipped_WAVECAR=gzipped_WAVECAR, backup=backup)
        else:
            my_wf = use_custodian(my_wf)

        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp)

        fw = self.lp.get_fw_by_id(fw_id=fw_id)
        with open(os.path.join(fw.launches[-1].launch_dir, "task_lobster.json")) as f:
            d = json.load(f)
        self._check_run(d, mode="lobsternormal")

        wf = self.lp.get_wf_by_fw_id(fw_id)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    def _single_Lobster_db_insertion(self, delete_wavecars=False, gzipped_output=False, gzipped_WAVECAR=False, fw_id=1,
                                     user_supplied_basis=None, backup=False, fake=True):
        # add the workflow

        structure = self.struct_si
        # with ScratchDir('.', copy_from_current_on_enter=True) as d:
        my_wf = get_wf_lobster(structure=structure, c={"vasp_cmd": VASP_CMD, "DB_FILE": DB_FILE},
                           user_kpoints_settings={"grid_density": 100}, delete_all_wavecars=delete_wavecars,
                           user_supplied_basis=user_supplied_basis, material_id="normal")
        if fake:
            my_wf = use_fake_vasp(my_wf, refs_dirs_si_vasp)
            my_wf = use_fake_lobster(my_wf, refs_dirs_si_lobster, gzipped_output=gzipped_output,
                                     gzipped_WAVECAR=gzipped_WAVECAR, backup=backup)
        else:
            my_wf = use_custodian(my_wf)

        self.lp.add_wf(my_wf)

        # run the workflow
        rapidfire(self.lp)
        d = self.get_task_collection(coll_name="lobster").find_one({"material_id": "normal"})
        self._check_run(d, mode="lobsternormal")

        wf = self.lp.get_wf_by_fw_id(fw_id)
        self.assertTrue(all([s == 'COMPLETED' for s in wf.fw_states.values()]))

    def test_single_vasp_lobster(self):
        # will always run fake vasp and fake lobster
        # These are mainly options for the fake lobster runs, they won't occur in the real workflows
        self._single_Vasp_Lobster(gzipped_output=True, gzipped_WAVECAR=True, fw_id=1)
        #self._single_Vasp_Lobster(gzipped_output=True, gzipped_WAVECAR=False, fw_id=3)
        #self._single_Vasp_Lobster(gzipped_output=False, gzipped_WAVECAR=False, fw_id=5)
        #self._single_Vasp_Lobster(gzipped_output=False, gzipped_WAVECAR=True, fw_id=7)
        #self._single_Vasp_Lobster(delete_wavecars=True, gzipped_output=False, gzipped_WAVECAR=False, fw_id=9)
        # include user supplied basis
        #self._single_Vasp_Lobster(delete_wavecars=False, gzipped_output=False, gzipped_WAVECAR=False, fw_id=11,
        #                          user_supplied_basis={"Si": "3s 3p"})
        #self._single_Vasp_Lobster(delete_wavecars=False, gzipped_output=False, gzipped_WAVECAR=False, fw_id=13,
        #                         backup=False)
        #self._single_Lobster_db_insertion(fake=True, fw_id=15)
        # integration test
        #if VASP_CMD and LOBSTER_CMD:
        #    self._single_Vasp_Lobster(fake=False, fw_id=17)

    #def tearDown(self):
    #    pass



class Test_get_all_possible_basis_combinations(PymatgenTest):
    def setUp(self):
        pass

    def test_getting_basis_combinations(self):
        # this basis is just for testing (not correct)
        min_basis = ['Li 1s 2s ', 'Na 1s 2s', 'Si 1s 2s']
        max_basis = ['Li 1s 2p 2s ', 'Na 1s 2p 2s', 'Si 1s 2s']
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        self.assertListEqual(combinations_basis,
                             [['Li 1s 2s', 'Na 1s 2s', 'Si 1s 2s'], ['Li 1s 2s', 'Na 1s 2s 2p', 'Si 1s 2s'],
                              ['Li 1s 2s 2p', 'Na 1s 2s', 'Si 1s 2s'], ['Li 1s 2s 2p', 'Na 1s 2s 2p', 'Si 1s 2s']])

        min_basis = ['Li 1s 2s']
        max_basis = ['Li 1s 2s 2p 3s']
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        self.assertListEqual(combinations_basis, [['Li 1s 2s'], ['Li 1s 2s 2p'], ['Li 1s 2s 3s'], ['Li 1s 2s 2p 3s']])

        min_basis = ['Li 1s 2s', 'Na 1s 2s']
        max_basis = ['Li 1s 2s 2p 3s', 'Na 1s 2s 2p 3s']
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        self.assertListEqual(combinations_basis,
                             [['Li 1s 2s', 'Na 1s 2s'], ['Li 1s 2s', 'Na 1s 2s 2p'], ['Li 1s 2s', 'Na 1s 2s 3s'],
                              ['Li 1s 2s', 'Na 1s 2s 2p 3s'], ['Li 1s 2s 2p', 'Na 1s 2s'],
                              ['Li 1s 2s 2p', 'Na 1s 2s 2p'], ['Li 1s 2s 2p', 'Na 1s 2s 3s'],
                              ['Li 1s 2s 2p', 'Na 1s 2s 2p 3s'], ['Li 1s 2s 3s', 'Na 1s 2s'],
                              ['Li 1s 2s 3s', 'Na 1s 2s 2p'], ['Li 1s 2s 3s', 'Na 1s 2s 3s'],
                              ['Li 1s 2s 3s', 'Na 1s 2s 2p 3s'], ['Li 1s 2s 2p 3s', 'Na 1s 2s'],
                              ['Li 1s 2s 2p 3s', 'Na 1s 2s 2p'], ['Li 1s 2s 2p 3s', 'Na 1s 2s 3s'],
                              ['Li 1s 2s 2p 3s', 'Na 1s 2s 2p 3s']])

        min_basis = ['Si 1s 2s 2p', 'Na 1s 2s']
        max_basis = ['Si 1s 2s 2p 3s', 'Na 1s 2s 2p 3s']
        combinations_basis = get_all_possible_basis_combinations(min_basis, max_basis)
        self.assertListEqual(combinations_basis, [['Si 1s 2s 2p', 'Na 1s 2s'], ['Si 1s 2s 2p', 'Na 1s 2s 2p'],
                                                  ['Si 1s 2s 2p', 'Na 1s 2s 3s'], ['Si 1s 2s 2p', 'Na 1s 2s 2p 3s'],
                                                  ['Si 1s 2s 2p 3s', 'Na 1s 2s'], ['Si 1s 2s 2p 3s', 'Na 1s 2s 2p'],
                                                  ['Si 1s 2s 2p 3s', 'Na 1s 2s 3s'],
                                                  ['Si 1s 2s 2p 3s', 'Na 1s 2s 2p 3s']])



if __name__ == "__main__":
    unittest.main()

