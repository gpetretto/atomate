import json
import json
import logging
import os
import shutil

from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import env_chk, get_meta_from_structure
from atomate.vasp.database import VaspCalcDb
from custodian import Custodian

from custodian.lobster.jobs import LobsterJob
from custodian.lobster.handlers import ChargeSpillingValidator, EnoughBandsValidator, LobsterFilesValidator, \
    CrashErrorHandler
from fireworks import FiretaskBase, explicit_serialize, FWAction
from fireworks.core.firework import Workflow
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from monty.json import jsanitize
from monty.shutil import compress_file
from pymatgen.core.structure import Structure
from pymatgen.io.lobster import Lobsterout, Lobsterin

LOBSTERINPUT_FILES = ["lobsterin"]
LOBSTEROUTPUT_FILES = ["lobsterout", "CHARGE.lobster", "COHPCAR.lobster", "COOPCAR.lobster", "DOSCAR.lobster",
                       "GROSSPOP.lobster", "ICOHPLIST.lobster", "ICOOPLIST.lobster", "lobster.out",
                       "projectionData.lobster"]

VASP_OUTPUT_FILES = ["OUTCAR", "vasprun.xml", "CHG", "CHGCAR", "CONTCAR", "INCAR", "KPOINTS", "POSCAR", "POTCAR",
                     "DOSCAR", "EIGENVAL", "IBZKPT", "OSZICAR", "PCDAT", "PROCAR", "REPORT", "WAVECAR", "XDATCAR"]

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
logger = logging.getLogger(__name__)


@explicit_serialize
class WriteLobsterinputfromIO(FiretaskBase):
    """
    will write lobsterin from POSCAR, INCAR, POTCAR
    Required Params:
        poscaraddress (str): address to POSCAR
        poscaraddress (str): address to POSCAR
        poscaraddress (str): address to POSCAR
        option (str): options as in Lobsterin.standard_calculations_from_vasp_files
    Optional Params:
        user_supplied_basis (dict): dictionary including the basis for each atom type
        user_lobsterin_settings (dict): dictionary that will be used to overwrite settings in Lobsterin dict
    """
    required_params = ["poscaraddress", "incaraddress", "potcaraddress", "option"]
    optional_params = ["user_supplied_basis", "user_lobsterin_settings"]

    def run_task(self, fw_spec):
        poscar = self.get("poscaraddress", "POSCAR")
        incar = self.get("incaraddress", "INCAR")
        potcar = self.get("potcaraddress", "POTCAR")
        option = self.get("option", "standard")
        user_supplied_basis = self.get("user_supplied_basis", None)
        if user_supplied_basis is None:
            lobsterinput = Lobsterin.standard_calculations_from_vasp_files(poscar, incar, potcar, option=option)
        else:
            lobsterinput = Lobsterin.standard_calculations_from_vasp_files(poscar, incar, None, option=option,
                                                                           dict_for_basis=user_supplied_basis)
        additional_input = self.get("user_lobsterin_settings", None)
        if additional_input:
            for key, parameter in additional_input.items():
                lobsterinput[key] = parameter

        lobsterinput.write_lobsterin("lobsterin")


# TODO: can we test this somehow?
@explicit_serialize
class RunLobster(FiretaskBase):
    """
    Starts the Lobster Job
    Optional params:
        lobster_cmd (str): command to run lobster
        gzip_output (bool): If true, output (except WAVECAR) will be gzipped.
        gzip_WAVECAR (bool): If true, WAVECAR will be gzipped
        strict_handlers_validators (bool): If True, the additional validator of the charge spilling will be used
    """

    required_params = []
    optional_params = ["lobster_cmd", "gzip_output", "gzip_WAVECAR", "strict_handlers_validators"]

    def run_task(self, fw_spec):
        lobster_cmd = env_chk(self.get('lobster_cmd'), fw_spec)
        gzip_output = self.get("gzip_output", False)
        gzip_WAVECAR = self.get("gzip_WAVECAR", False)
        if gzip_WAVECAR:
            # TODO: include the changes in the lobster job correctly
            add_files_to_gzip = VASP_OUTPUT_FILES
        else:
            add_files_to_gzip = VASP_OUTPUT_FILES.remove("WAVECAR")

        strict_handlers_validators = self.get("strict_handlers_validators", False)
        # LobsterJob gzips output files, Custodian would gzip all output files (even slurm)
        jobs = [LobsterJob(lobster_cmd=lobster_cmd, output_file="lobster.out", stderr_file="std_err_lobster.txt",
                           gzipped=gzip_output, add_files_to_gzip=add_files_to_gzip)]
        handlers = [CrashErrorHandler(output_filename="lobsterout")]
        # TODO: should ChargeSpillingValdator be included? This might lead to confusion!
        if strict_handlers_validators:
            validators = [ChargeSpillingValidator(output_filename="lobsterout"), LobsterFilesValidator(),
                          EnoughBandsValidator(output_filename="lobsterout")]
        else:
            validators = [LobsterFilesValidator(),
                          EnoughBandsValidator(output_filename="lobsterout")]
        c = Custodian(handlers=handlers, jobs=jobs, validators=validators, gzipped_output=False, max_errors=5)
        c.run()


@explicit_serialize
class LobsterRunToDb(FiretaskBase):
    """
    Adds Lobster Calculation to Database.Uses current directory unless you
    specify calc_dir or calc_loc.
    Optional params:
        calc_dir (str): path to dir (on current filesystem) that contains VASP
            output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        additional_fields (dict): dict of additional fields to add
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
    """
    # TODO: which ones are required, which are optional
    optional_params = ["calc_dir", "calc_loc", "additional_fields", "db_file", "fw_spec_field"]

    def run_task(self, fw_spec):

        vasp_calc_dir = self.get("calc_dir", None)
        # TODO: ist this even correct?
        vasp_calc_loc = get_calc_loc(self["calc_loc"], fw_spec["calc_locs"]) if self.get("calc_loc") else {}

        # get the directory that contains the Lobster dir to parse
        calc_dir = os.getcwd()
        # parse the Lobster directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))
        if os.path.exists("POSCAR.gz"):
            struct = Structure.from_file("POSCAR.gz")
        elif os.path.exists("POSCAR"):
            struct = Structure.from_file("POSCAR")
        else:
            raise ValueError("POSCAR.gz/POSCAR does not exist")

        # use lobsterout from pymatgen to get all the information
        if os.path.exists("lobsterout.gz"):
            Lobsterout_here = Lobsterout("lobsterout.gz")
        elif os.path.exists("lobsterout"):
            Lobsterout_here = Lobsterout("lobsterout")
        else:
            raise ValueError("lobsterout.gz/lobsterout does not exist")
        # will get document from Lobsterout
        task_doc = {}
        task_doc["output"] = jsanitize(Lobsterout_here.get_doc())
        if os.path.exists("lobsterin.gz"):
            task_doc["input"] = jsanitize(Lobsterin.from_file("lobsterin.gz"))
        elif os.path.exists("lobsterin"):
            task_doc["input"] = jsanitize(Lobsterin.from_file("lobsterin"))
        else:
            raise ValueError("lobsterin.gz/lobsterin does not exist")
        if os.path.exists("lobsterin.orig.gz"):
            task_doc["orig_input"] = jsanitize(Lobsterin.from_file("lobsterin.orig.gz"))
        elif os.path.exists("lobsterin.orig"):
            task_doc["orig_input"] = jsanitize(Lobsterin.from_file("lobsterin.orig"))
        # save custodian details
        try:
            with open("custodian.json", "r") as f:
                custodian_details = json.load(f)
            task_doc["custodian"] = jsanitize(custodian_details)
        except:
            pass
        additional_fields = self.get("additional_fields", {})
        if additional_fields:
            task_doc.update(additional_fields["additional_fields"])

        task_doc.update(get_meta_from_structure(struct))
        if vasp_calc_dir != None:
            task_doc["vasp_dir_name"] = vasp_calc_dir
        else: 
            task_doc["vasp_dir_name"] = vasp_calc_loc["path"]
        task_doc["dir_name"] = os.getcwd()

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update(fw_spec[self.get("fw_spec_field")])

        task_doc["state"] = 'successful'

        task_doc = jsanitize(task_doc)
        # get the database connection
        db_file = env_chk(self.get('db_file'), fw_spec)

        # db insertion or taskdoc dump
        if not db_file:
            with open("task_lobster.json", "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            db = VaspCalcDb.from_db_file(db_file, admin=True)
            db.collection = db.db["lobster"]
            db.collection.insert(task_doc)
            logger.info("Lobster calculation is complete.")
        return FWAction()


#TODO: adapat it to make it usable with more complex workflows (more than 1 lobster calc)
def use_fake_lobster(original_wf, ref_dirs, params_to_check=None, gzipped_output=False, gzipped_WAVECAR=False,
                     backup=False):
    """
    NOT TESTED YET
    Replaces all tasks with "RunLobster" to be RunLobsterFake. Thus, we do not
    actually run Lobster but copy pre-determined inputs and outputs.
    Args:
        original_wf (Workflow)
        ref_dirs (dict): key=firework name, value=path to the reference lobster calculation directory
        params_to_check (list): optional list of lobsterin parameters that are checked
        gzipped_output (bool): if true, everything except WAVECAR will be gzipped
        gzipped_WAVECAR (bool): if true, WAVECAR will be gzipped
        backup (bool): if true, lobsterin will be saved to lobsterin.orig
    Returns:
        Workflow
    """
    if not params_to_check:
        params_to_check = ["basisSet", "cohpGenerator"]
    for idx_fw, fw in enumerate(original_wf.fws):
        for job_type in ref_dirs.keys():
            if job_type in fw.name:
                for idx_t, t in enumerate(fw.tasks):
                    print(t)
                    print(job_type)
                    if "RunLobster" in str(t):
                        original_wf.fws[idx_fw].tasks[idx_t] = \
                            RunLobsterFake(ref_dir=ref_dirs[job_type],
                                           params_to_check=params_to_check, gzipped_output=gzipped_output,
                                           gzipped_WAVECAR=gzipped_WAVECAR, backup=backup)

    return original_wf


@explicit_serialize
class RunLobsterFake(FiretaskBase):
    """
     Lobster Emulator; not yet tested
     Required params:
         ref_dir (string): Path to reference lobster run directory with input files in the folder
            named 'inputs' and output files in the folder named 'outputs'.
     Optional params:
         params_to_check (list): optional list of lobsterin parameters to check
         check_lobsterin (bool): whether to confirm the lobsterin params (default: True)
     """
    required_params = ["ref_dir"]
    optional_params = ["params_to_check", "check_lobsterin", "gzipped_output", "gzipped_WAVECAR", "backup"]

    def run_task(self, fw_spec):
        self._verify_inputs()
        self._clear_inputs()
        self._generate_outputs()

    def _verify_inputs(self):
        user_lobsterin = Lobsterin.from_file(os.path.join(os.getcwd(), "lobsterin"))

        # Carry out some BASIC tests.

        # Check lobsterin
        if self.get("check_lobsterin", True):
            #TODO understand this class better
            ref_lobsterin = Lobsterin.from_file(os.path.join(self["ref_dir"], "inputs", "lobsterin"))
            params_to_check = self.get("params_to_check", [])
            defaults = {"basisSet": "pbeVaspFit2015", "cohpEndEnergy": 5.0}
            for p in params_to_check:
                if user_lobsterin.get(p, defaults.get(p)) != ref_lobsterin.get(p, defaults.get(p)):
                    raise ValueError("lobsterin value of {} is inconsistent!".format(p))

        logger.info("RunLobsterFake: verified inputs successfully")

    def _clear_inputs(self):
        for x in ["lobsterin"]:
            p = os.path.join(os.getcwd(), x)
            if os.path.exists(p):
                os.remove(p)

    def _generate_outputs(self):
        # generate lobsterin input
        # pretend to have run lobster by copying pre-generated outputs from reference dir to cur dir
        output_dir = os.path.join(self["ref_dir"], "outputs")
        for file_name in os.listdir(output_dir):
            full_file_name = os.path.join(output_dir, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, os.getcwd())

        if self.get("gzipped_WAVECAR", False):
            for file in LOBSTEROUTPUT_FILES:
                if os.path.exists(file):
                    compress_file(file, compression="gz")
            for file in LOBSTERINPUT_FILES:
                if os.path.exists(file):
                    compress_file(file, compression="gz")
            if self.get("backup", False):
                if os.path.exists("lobsterin.orig"):
                    compress_file("lobsterin.orig", compression="gz")
            for file in VASP_OUTPUT_FILES:
                if self.get("gzipped_WAVECAR", False):
                    if os.path.exists(file):
                        compress_file(file, compression="gz")
                else:
                    if os.path.exists(file) and file != 'WAVECAR':
                        compress_file(file, compression="gz")

        logger.info("RunLobsterFake: ran fake lobster, generated outputs")
