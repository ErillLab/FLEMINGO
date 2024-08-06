from . import utils
from . import config
from .genome_object import Genome
from Bio.Seq import Seq
from Markov_DNA import MCM
import math
from .exporters.exporter import Exporter
from .organism_factory import OrganismFactory
import numpy as np

class Scanner:

    def __init__(self) -> None:
        """ Initializes the scanner class, which contains all the methods for 
           performing genomic scans with composite motifs.
        """
        self.placements = []
        self.threshold_scores = []
        self.threshold = 0
        self.M = 1.5

    def set_up(self, cfg_file: str = config.JSON_CONFIG_FILENAME) -> None:
        """ Reads in the config file and establishes global variables
        """
        config.set_up(cfg_file=cfg_file)

    def load(self) -> None:
        """ Initializes the scanner:
            - Loads up genome
            - Loads up model
            - Computes window size for scanning
            - Determines score threshold
        """
        # genome and model loading
        self.genome = Genome(config.INPUT_GENOME_FILE, format=config.INPUT_TYPE)
        self.model = OrganismFactory.import_scanner_model(config.INPUT_MODEL_FILE)
        
        if config.VERBOSE:
            self.model.print()

        # determine window size, using the method from the chain model to estimate its
        # placement length within a confidence interval
        self.WS = self.model.get_max_length(config.LENGTH_CI)

        # determine threshold by generating N random sequences of length equal to window
        # and evaluating the model placement in them, then using the significance value 
        # to determine the threshold
        self.threshold = self._get_threshold_score()

    def scan(self) -> None:
        """ Implements the main scan process. Scans the specified regions 
            by the user (ALL THE GENOME, INTERGENIC REGIONS or REGIONS OF
            DETERMINED SIZE UPSTREAM THE GENES) searching for significant 
            placements.
                - INTERGENIC: intergenic regions.
                - UPSTREAM_REGIONS: regions of determined size upstream 
                  the genes.
                - WHOLE: scans all the genome.

            The scanning is performed using the sliding window method. The 
            size of the window is determined based on the provided composite 
            motif model and the step size is half of the window size. In 
            each of the windows the optimal placement is computed, and in 
            case that outperforms the threshold it is considered as a significant 
            placement. For significant placements the nearest genes are obtained.
        """
        # Determine the regions to be scanned. Region defined as [START, END]
        # whether we only scan intergenic regions (INTERGENIC)
        # scan gene upstream regions (UPSTREAM_REGIONS) [both of these apply to GenBank files]
        # or scan ALL the genome [typically used for FASTA files]
        if config.SCAN_MODE == "INTERGENIC":
            regions = self.genome.inter_regions        #list of regions
        
        elif config.SCAN_MODE == "UPSTREAM_REGIONS":
            regions = self.genome.upstream_regions     #list of regions
        elif config.SCAN_MODE == "WHOLE":
            regions = [[0, self.genome.length]]        #only one region
        else:
            print("ERROR: Invalid scan method selected.")
            print("Select between WHOLE, INTERGENIC or UPSTREAM_REGIONS modes")
            exit(1)

        # variables to prevent adding multiple times same placement
        # across overlapping windows
        prev_score = -1
        prev_loc = -1

        # for each region to scan
        for reg in regions:
            # determine length of current region
            # ex: intergenic region has 700 bp
            size = reg[1]-reg[0]

            # obtain total size to be scanned, as the 
            # number of windows to be used (size // self.WS)+1)
            # times the length of the window
            # ex: the 700 bp will be encased in a segment made by a multiple of WS
            # ex: for WS = 500; 700 // 500 = 1 [1.4] + 1 => 2 --> 2*500 = 1000
            # ex: we would be scanning 1000 bp
            aux_N = ((size // self.WS)+1)*self.WS

            # start position for scan (offset) is obtained as the
            # center of the region (reg[0] + size//2), minus half
            # the length of the scanned region, ensuring that the
            # scanned segment encases the window in the middle
            # ex: assume region starts at genomic position 900
            # ex: 900 + 700//2 = 900 + 350 = 1250
            # ex: 1250 - 1000//2 = 750 <-- our offset starts 150 (300//2) upstream of region
            offset = reg[0] + size//2 - aux_N//2

            # out of bounds checks
            if offset < 0:    
                offset = 0
            if aux_N > self.genome.length:
                aux_N = self.genome.length

            # scanning loop, using WS fragments and a step size WS // 2
            for n in range(0, aux_N - (self.WS // 2), self.WS // 2):
                window = str(self.genome.seq[n+offset:n+offset+self.WS])
                if len(window) < self.WS // 2:
                    break
                # get optimal placement in both strands (if required)
                for_placement = self.model.get_placement(window)
                if config.SCAN_REVERSE:
                    window = str(Seq(window).reverse_complement())
                    rev_placement = self.model.get_placement(window)

                # pick optimal placement (for both strands)
                if not config.SCAN_REVERSE or for_placement.energy >= rev_placement.energy:
                    placement = for_placement
                    placement.strand = "+"
                else:
                    placement = rev_placement
                    placement.strand = "-"
                    # update coordinates for placement on the reverse strand
                    placement = utils.reverse_placement(placement, len(window))

                # if placement does not make the cut or it is the same placement as in previous window
                if placement.energy < self.threshold or (prev_score == placement.energy and prev_loc == offset + n + placement.start_pos()):
                    continue

                # if placement is above the threshold
                # identify nearest genes and associate to placement
                l_gene, r_gene = self.genome.get_nearest_genes(offset + n + placement.start_pos(), offset + n + placement.end_pos())
                placement.r_gene = r_gene
                placement.l_gene = l_gene

                # define window coordinates that generated the placement
                placement.window=[n+offset, n+offset+self.WS]    

                if config.VERBOSE:
                    utils.print_placement_info(placement)    

                # obtain genomic location for placement    
                placement.add_offset(n + offset)

                # compute the p-value for the placement
                placement.p_value = self._compute_p_value(rev_placement.energy)
                self.placements.append(placement)


    def export(self, filename: str, format=Exporter.E_BED) -> None:
        """ Exports the results (significant placements information) 
            to an output file with the specified file format. The 
            posible formats are:

                - JSON
                - BED
                - CSV
                - GFF3
        """
        exporter = Exporter(filename, format=format)
        if format == Exporter.E_JSON:
            exporter.export(self.genome, self.placements, self.model)
        else:
            exporter.export(self.genome, self.placements)      

    def export_genes(self, filename) -> None:
        """ Exports the information of the nearest genes from
            each significant placement to a CSV file
        """
        genes = []
        # gets the nearest genes from each placement and appends them to a list
        for placement in self.placements:
            # gets the information of the left gene
            gene = placement.l_gene
            if gene != None:
                if gene.score == None:
                    genes.append(gene)
                if gene.score == None or gene.score < placement.energy:
                    gene.score = placement.energy
                    gene.p_value = placement.p_value

            # gets the information of the left gene
            gene = placement.r_gene
            if gene != None:
                if gene.score == None:
                    genes.append(gene)
                if gene.score == None or gene.score < placement.energy:
                    gene.score = placement.energy
                    gene.p_value = placement.p_value

        # writes CSV file with genes information
        with open(filename, 'w') as f:
            print("Tag", "Name", "Start", "End", "Strand", "Score", "p-value", "Protein_Accession", "Product", sep=",", file=f)
            for gene in genes:
                print(gene.tag, gene.name, gene.location.start, gene.location.end, gene.strand, gene.score, gene.p_value, gene.protein_id, '"'+str(gene.product)+'"', sep=",", file=f)
    
    def refine(self) -> None:
        """ Performs iterative refinement of the p-value estimated for a placement
            using a sampling size of SAMPLE_SIZE (number of pseudosequences) a terminal
            p-value to be reached (THRESHOLD_P_VALUE) and a Markov Model order (k).
        """
        # estimate the best MM order (k) to used based on user preference
        # and length of the sequence it is trained on [shorter sequences lead to
        # lower order Markov Models]
        aux_k = self.infer_k(self.placements[0].dna_sequence)
        # selects the minimun k between the specified one by the user or the inferred one
        k = min(config.REFINEMENT_K, aux_k)
        if config.VERBOSE:
            print("Refinement k", k)
        
        N = 10    # initial sampling size
        i=1       # iteration
        while N <= config.SAMPLE_SIZE:      # while not reaching terminal sampling size
            if config.VERBOSE:
                print("Refinement iteration", i)
            self._refinement_it(N, k)       # make refinement iteration
            if config.VERBOSE:
                print("Remaining placements:", len(self.placements))
            N*=10
            i+=1

    def _refinement_it(self, N: int, k: int) -> None:
        """ Performs a p-value refinement iteration
            using a MM of order k, and a sampling size of N

            In each refinement iteration it is computed a new
            score threshold for each significant placement. This 
            new threshold is computed using N pseudo-replicates
            generated by MM of order k. 
        """
        # number of pseudoreplicates with score > observed to reach p-value
        threshold = int(config.THRESHOLD_P_VALUE * N)

        new_placements = []

        # for each window with a significant placement
        for placement in self.placements:
            # obtain the sequence (window) that generated the placement
            seq = placement.dna_sequence
            # generate N pseudo-replicates of the window sequence using MM of order k
            mcm = MCM(k)
            mcm.train(seq)
            seqs = mcm.generate(N)
            scores = []
            # perform placements on all pseudo-replicates
            for seq in seqs:
                aux_placement = self.model.get_placement(seq)
                scores.append(aux_placement.energy)
            # sort scores and determine if placement score is above threshold
            scores.sort()
            score = scores[-threshold]
            if placement.energy >= score:
                new_placements.append(placement)

        # update list of placements with the ones making the cut
        self.placements = new_placements

    def _get_threshold_score(self) -> float:
        """ Computes initial score threshold for a given model.
            Estimates placement scores on SAMPLE_SIZE random sequences 
            and computes p-value from this.
        """
        score = 0
        # compute the position on the scrores array for initial score threshold 
        threshold = int(config.SAMPLE_SIZE*config.THRESHOLD_P_VALUE)
        scores = []
        seqs = []

        # generate SAMPLE_SIZE random sequences
        for _ in range(config.SAMPLE_SIZE):
            new_seq = np.random.choice(["A","C","G","T"], p=[0.25, 0.25, 0.25, 0.25], size=self.WS)
            new_seq = "".join(new_seq)
            seqs.append(new_seq)

        # for each sequence, obtain optimal placement
        for seq in seqs:
            placement = self.model.get_placement(seq)
            scores.append(placement.energy)

        # sort scores and grab the score that corresponds to p-value limit
        scores.sort()
        self.threshold_scores = scores
        score = scores[-threshold]
        return score

    def infer_k(self, seq: str) -> int:
        """ Computes Markov Model order based on sequence length 

            For a given value of k, there are 4**k k-mers. We seek a value of k 
            such that the expected value of the number of occurrences of a k-mer 
            is at least M. Under the assumption that all the 4**k k-mers are equally 
            likely, the expected value of the count (number of occurrences) of any 
            k-mer is L/(4**k), where L is the length of the sequence on which we 
            train the Markov Model. Hence, we require:

            L/(4**k) >= M --> (4**k) <= L/M --> k <= log4(L/M)
        """
        L = len(seq)
        return int(math.log(L/self.M, 4))-1
    
    def _compute_p_value(self, score: float) -> float:
        """ Computes the p-value for a given score based on the array
            of scores generated in the threshold determintation.

            Determines the p-value for a placement that has passed the
            threshold.
        """
        # get the number of scores above the placement score
        count = len(list(filter(lambda x: x >= score, self.threshold_scores)))
        # compute the p-value
        return count/len(self.threshold_scores)
