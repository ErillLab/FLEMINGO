from objects.scanner import Scanner
import argparse
import objects.config as config
from Markov_DNA import MCM

# Main program for command line call
if __name__ == "__main__":

    # Parse command line arguments (only config file name expected)
    parser = argparse.ArgumentParser(
                            prog='scanner.py',
                            description='Runs the FLEMINGO scanner to detect instances of composite motifs in genome sequences',
                            epilog='Specify the config file name (default: python main.py -c scanner_config.json')
    parser.add_argument("-c", "--config", type=str, default="scanner_config.json")
    args = parser.parse_args()
    print(args.config)

    # Instantiates the scanner object, which handles the scanning process
    scanner = Scanner()
    # Reads in the configuration file and sets global variables
    scanner.set_up(args.config)
    # Loads the genome and composite motif model with file names specified in configuration file
    scanner.load()

    # Runs the scanner
    scanner.scan()
    # Reports the initial score threshold for significant composite motif placements
    # inferred from placements on randomize sequences
    print(scanner.threshold)

    # Exports unrefined (no_ref) scan results, appending relevant suffixes
    # Export to GFF3 format (all placements), 
    scanner.export_GFF3(config.OUTPUT_FILE.replace(".json", "_no_ref.gff3"))
    # Export list of genes near placements (putatively regulated genes)
    scanner.export_genes(config.OUTPUT_FILE.replace(".json", "_no_ref_genes.csv"))
    # Export list of all significant placements (score, sequence...) and information on putatively regulated genes
    scanner.export_CSV(config.OUTPUT_FILE.replace(".json", "_no_ref.csv"))
    # Iterative refinement of p-value associated with each placement
    # For each significant placement, generate N pseudoreplicates with k-order Markov Model to estimate p-value
    scanner.refine(p_value=config.THRESHOLD_P_VALUE, n=1000, k=7)    # <-- make this dependent on user input (n=1000, k=7 )
    # Export scan results after refinement, appending relevant suffixes
    scanner.export_JSON(config.OUTPUT_FILE)
    scanner.export_BED(config.OUTPUT_FILE.replace(".json", ".bed"))
    scanner.export_genes(config.OUTPUT_FILE.replace(".json", "_genes.csv"))
    scanner.export_GFF3(config.OUTPUT_FILE.replace(".json", ".gff3"))
    scanner.export_CSV(config.OUTPUT_FILE.replace(".json", ".csv"))
