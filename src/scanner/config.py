import time
import shutil
import FMS.utils as utils


JSON_CONFIG_FILENAME = "config.json"

def set_up(cfg_file=JSON_CONFIG_FILENAME):
    """Reads configuration file and sets up all program variables

    """

    # specify as global variable so it can be accesed in local
    # contexts outside setUp
    global RUN_MODE
    global GEN_INFO_TO_STDOUT
    global INPUT_GENOME_FILE
    global INPUT_ORGANISM_FILE
    global LENGTH_CI
    global OUTPUT_FILE
    global INPUT_TYPE
    global VERBOSE
    global SCAN_MODE
    global SCAN_REVERSE
    global THRESHOLD_P_VALUE
    global UPSTREAM_REGION_SIZE
    global SAMPLE_SIZE

    # MPI variables
    global comm
    global rank
    global p

    config = utils.read_json_file(cfg_file)

    # Check the config settings and set up MPI if running in parallel
    RUN_MODE = config["main"]["RUN_MODE"]
    if RUN_MODE == "parallel":
        from mpi4py import MPI  # mpi4py is only imported if needed
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        p = comm.Get_size()
        # Make all processes wait for process 0 to complete the check
        comm.Barrier()
    elif RUN_MODE == "serial":
        comm, rank, p = None, None, None
    else:
        raise ValueError('RUN_MODE should be "serial" or "parallel".')

    if utils.i_am_main_process():
        print('\n==================\nRUN_MODE: {}\n==================\n'.format(RUN_MODE))

    GEN_INFO_TO_STDOUT = config["main"]["GEN_INFO_TO_STDOUT"]
    INPUT_GENOME_FILE = config["main"]["INPUT_GENOME_FILE"]
    INPUT_ORGANISM_FILE = config["main"]["INPUT_ORGANISM_FILE"]
    LENGTH_CI = config["main"]["LENGTH_CI"]
    OUTPUT_FILE = config["main"]["OUTPUT_FILE"]
    INPUT_TYPE = utils.get_file_type(INPUT_GENOME_FILE)
    VERBOSE = config["main"]["VERBOSE"]
    SCAN_MODE = config["main"]["SCAN_MODE"] # "INTERGENIC", "PROMOTER", "WHOLE"
    if SCAN_MODE == "UPSTREAM_REGIONS":
        UPSTREAM_REGION_SIZE = config["main"]["UPSTREAM_REGION_SIZE"]
    SCAN_REVERSE = config["main"]["SCAN_REVERSE"]
    THRESHOLD_P_VALUE = config["main"]["THRESHOLD_P_VALUE"]
    SAMPLE_SIZE = config["main"]["SAMPLE_SIZE"]

    # Throw config on a file
    if utils.i_am_main_process():
        print_config_json(config["main"], "Main Config")

def print_config_json(config: dict, name: str) -> None:
    """Print the config file on std out and send it to a file.
    It is useful so we can know which was the configuration on every run

    Args:
        config: Configuration file to print
        name: Title for the configuration file
        path: File to export the configuration info
    """
    utils.print_ln("{}:".format(name))

    for key in config.keys():
        utils.print_ln("{}: {}".format(key, config[key]))
    utils.print_ln("\n")

