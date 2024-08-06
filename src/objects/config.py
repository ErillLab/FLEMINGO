from . import utils

JSON_CONFIG_FILENAME = "scanner_config.json"

def set_up(cfg_file=JSON_CONFIG_FILENAME):
    """ Reads configuration file and sets up all program variables

    """
    # specify as global variable so it can be accessed in local
    # contexts outside setUp
    
    # filename of the file which contains the genome
    global INPUT_GENOME_FILE
    # filename of the file which contains the model
    global INPUT_MODEL_FILE
    # confidence interval used for window size computation
    global LENGTH_CI
    # output filename
    global OUTPUT_FILE
    # input genome format (genbank or fasta)
    global INPUT_TYPE   
    # output some information 
    global VERBOSE
    # determines the scan mode (WHOLE, INTERGENIC, UPSTREAM_REGIONS)
    global SCAN_MODE
    # in case of scanning upstream regions of the genes, the length of those regions
    global UPSTREAM_REGION_SIZE
    # determines if the reverse strand is scanned
    global SCAN_REVERSE
    # p-value used for score threshold computations
    global THRESHOLD_P_VALUE
    # number of sequences used as sample to compute score threshold
    global SAMPLE_SIZE
    # perform results refinement
    global REFINEMENT
    # Markov Model order for pseudo-replicates generation in refinement phase
    global REFINEMENT_K
    # Output file formats
    global OUTPUT_FORMATS

    config = utils.read_json_file(cfg_file)

    # read the parameters
    INPUT_GENOME_FILE = config["main"]["INPUT_GENOME_FILE"]
    INPUT_MODEL_FILE = config["main"]["INPUT_MODEL_FILE"]
    LENGTH_CI = config["main"]["LENGTH_CI"]
    OUTPUT_FILE = config["main"]["OUTPUT_FILE"]
    INPUT_TYPE = utils.get_file_type(INPUT_GENOME_FILE)
    VERBOSE = config["main"]["VERBOSE"]
    SCAN_MODE = config["main"]["SCAN_MODE"] # "INTERGENIC", "UPSTREAM_REGIONS", "WHOLE"
    if SCAN_MODE == "UPSTREAM_REGIONS":
        UPSTREAM_REGION_SIZE = config["main"]["UPSTREAM_REGION_SIZE"]
    SCAN_REVERSE = config["main"]["SCAN_REVERSE"]
    THRESHOLD_P_VALUE = config["main"]["THRESHOLD_P_VALUE"]
    SAMPLE_SIZE = config["main"]["SAMPLE_SIZE"]
    REFINEMENT = config["main"]["REFINEMENT"]
    if REFINEMENT:
        REFINEMENT_K = config["main"]["REFINEMENT_K"]
    OUTPUT_FORMATS = config["main"]["OUTPUT_FORMATS"]
