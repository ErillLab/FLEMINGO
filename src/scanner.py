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

    # # Exports unrefined (no_ref) scan results, appending relevant suffixes
    # # Export to GFF3 format (all placements), 
    # scanner.export(config.OUTPUT_FILE.replace(".json", "_no_ref.gff3"), format="JSON")
    # # Export list of genes near placements (putatively regulated genes)
    # scanner.export_genes(config.OUTPUT_FILE.replace(".json", "_no_ref_genes.csv"))
    # # Export list of all significant placements (score, sequence...) and information on putatively regulated genes
    # scanner.export(config.OUTPUT_FILE.replace(".json", "_no_ref.csv"), format="JSON")

    # Iterative refinement of p-value associated with each placement
    # For each significant placement, generate N pseudoreplicates with k-order Markov Model to estimate p-value
    if config.REFINEMENT:
        scanner.refine()

    
    # Export scan results, appending relevant suffixes
    for format in config.OUTPUT_FORMATS:
        scanner.export(config.OUTPUT_FILE + "." + format.lower(), format=format)
    # scanner.export(config.OUTPUT_FILE + ".bed", format="BED")
    # scanner.export(config.OUTPUT_FILE + ".gff3", format="GFF3")
    # scanner.export(config.OUTPUT_FILE + ".csv", format="CSV")

    # Export gene list
    scanner.export_genes(config.OUTPUT_FILE + "_genes.csv")
