from objects.scanner import Scanner
import argparse
import FMS.config as config
from Bio import SeqIO
from Markov_DNA import MCM

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
                            prog='main.py',
                            description='What the program does',
                            epilog='Text at the bottom of help')
    parser.add_argument("-c", "--config", type=str, default="config.json")
    args = parser.parse_args()
    print(args.config)
    scanner = Scanner()
    scanner.set_up(args.config)
    scanner.load()

    scanner.scan()
    print(scanner.threshold)
    scanner.export_GFF3(config.OUTPUT_FILE.replace(".json", "_no_ref.gff3"))
    scanner.export_genes(config.OUTPUT_FILE.replace(".json", "_no_ref_genes.csv"))
    scanner.export_CSV(config.OUTPUT_FILE.replace(".json", "_no_ref.csv"))
    # scanner.refine(p_value=0.1, n=1000, k=7)
    scanner.refine(p_value=config.THRESHOLD_P_VALUE, n=1000, k=7)
    scanner.export()
    scanner.export_BED(config.OUTPUT_FILE.replace(".json", ".bed"))
    scanner.export_genes(config.OUTPUT_FILE.replace(".json", "_genes.csv"))
    scanner.export_GFF3(config.OUTPUT_FILE.replace(".json", ".gff3"))
    scanner.export_CSV(config.OUTPUT_FILE.replace(".json", ".csv"))
