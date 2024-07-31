import FMS.utils as utils
import FMS.exporter as exporter
import FMS.config as config
from FMS.objects.genome_object import Genome
from FMS.objects.chain_model_object import ChainModel
from FMS.objects.generator import Generator
from Bio.Seq import Seq
from Markov_DNA import MCM
import math
from FMS.objects.exporters.exporter import Exporter

class Scanner:

    def __init__(self, threshold_p_value=0.005, sample_size=10000):
        """Initalizes the scanner class, which contains all the methods for 
           performing genomic scans with composite motifs.
        """
        # self.placements_forw = []
        # self.placements_back = []
        self.placements = []
        self.threshold_scores = []
        self.threshold = 0
        self.threshold_p_value = threshold_p_value
        self.sample_size = sample_size
        self.NUCLEOTIDES = "ACGT"
        self.M = 1.5

    def set_up(self, cfg_file):
        """ Reads in the config file and establishes global variables
        """
        config.set_up(cfg_file=cfg_file)

    def load(self):
        """ Initalizes the scanner:
            - Loads up genome
            - Loads up model
            - Computes window size for scanning
            - Determines score threshold
        """

        # genome and model loading
        self.genome = Genome(config.INPUT_GENOME_FILE, format=config.INPUT_TYPE)
        self.model = ChainModel.import_model(config.INPUT_ORGANISM_FILE)
        
        if config.VERBOSE:
            self.model.print()

        # determine window size, using the method from the chain model to estimate its
        # placement length within a confidence interval
        self.WS = self.model.get_max_length(config.LENGTH_CI)
        print(self.WS)

        # determine threshold by generating N random sequences of length equal to window and
        # and evaluating the model placement in them, then using the significance value to determine the
        # threshold
        self.threshold_p_value = config.THRESHOLD_P_VALUE
        self.sample_size = config.SAMPLE_SIZE
        self.threshold = self.get_threshold_score(self.WS, N=self.sample_size, p_value=self.threshold_p_value)

    def scan(self):
        """ Implements the main scan process
        """

        # Determine the regions to be scanned. Region defined as [START, END]
        # whether we only scan intergenic regions (INTERGENIC)
        # scan gene upstream regions (UPSTREAM_REGIONS) [both of these apply to GenBank files]
        # or scan ALL the genome [typically used for FASTA files]
        if config.SCAN_MODE == "INTERGENIC":
            regions = self.genome.inter_regions        #list of regions
        
        elif config.SCAN_MODE == "UPSTREAM_REGIONS":
            regions = self.genome.upstream_regions     #list of regions
        else:
            regions = [[0, self.genome.length]]        #only one region

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

            #out of bounds checks
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

                #pick optimal placement (for both strands)
                if not config.SCAN_REVERSE or for_placement.energy >= rev_placement.energy:
                    placement = for_placement
                    placement.strand = "+"
                else:
                    placement = rev_placement
                    placement.strand = "-"
                    # update coordinates for placement on the reverse strand
                    placement = utils.reverse_placement(placement, len(window))

                # if placement does not make the cut or it is the same placement as in previous window
                if placement.energy < self.threshold or (prev_score == placement.energy and prev_loc == offset + n + for_placement.recognizers_positions[0][0]):
                    continue
                        

                # if placement is above the threshold
                # identify nearest genes and associate to placement
                l_gene, r_gene = self.genome.get_nearest_genes(offset + n + placement.recognizers_positions[0][0], offset + n + placement.recognizers_positions[-1][1], strand=placement.strand)
                placement.r_gene = r_gene
                placement.l_gene = l_gene
                if config.VERBOSE:
                    utils.print_placement_info(placement, l_gene, r_gene, n + offset)
                
                # obtain genomic location for placement    
                placement.add_offset(n + offset)
                # define window coordinates that generated the placement
                placement.window=[n+offset, n+offset+self.WS]
                # compute the p-value for the placement
                placement.p_value = self.compute_p_value(rev_placement.energy)
                self.placements.append(placement)
        return 0
    
    def export2(self):
        # for placement in self.placements_forw:
        #     exporter.export_placement(placement.window[0], placement.window[1], placement, None, None, strand="+")

        # for placement in self.placements_back:
        #     exporter.export_placement(placement.window[0], placement.window[1], placement, None, None, strand="-")

        for placement in self.placements_back:
            exporter.export_placement(placement.window[0], placement.window[1], placement, None, None, strand=placement.strand)

        exporter.export_results(self.genome.length, self.model, self.WS, self.genome.inter_regions)

    # def export_BED(self, filename):
    #     bed_exporter = Exporter(filename, format=Exporter.E_BED)
    #     bed_exporter.export(self.genome, self.placements_forw, self.placements_back)

    # def export_GFF3(self, filename):
    #     gff3_exporter = Exporter(filename, format=Exporter.E_GFF3)
    #     gff3_exporter.export(self.genome, self.placements_forw, self.placements_back)

    # def export_CSV(self, filename):
    #     csv_exporter = Exporter(filename, format=Exporter.E_CSV)
    #     csv_exporter.export(self.genome, self.placements_forw, self.placements_back)

    def export(self, filename, format=Exporter.E_BED):
        fexporter = Exporter(filename, format=format)
        fexporter.export(self.genome, self.placements, [])
        # fexporter.export(self.genome, self.placements_forw, self.placements_back)
        

    def export_genes(self, filename):
        genes = []
        for placement in self.placements:
            gene = placement.l_gene
            if gene != None and not gene in genes:
                gene.score = placement.energy
                gene.p_value = placement.p_value
                genes.append(gene)
            gene = placement.r_gene
            if gene != None and not gene in genes:
                gene.score = placement.energy
                gene.p_value = placement.p_value
                genes.append(gene)
        # for placement in self.placements_back:
        #     gene = placement.l_gene
        #     if gene != None and not gene in genes:
        #         genes.append(gene)
        #     gene = placement.r_gene
        #     if gene != None and not gene in genes:
        #         genes.append(gene)

        with open(filename, 'w') as f:
            print("Tag", "Name", "Start", "End", "Strand", "Score", "p-value", "Protein_Accession", "Product", sep=",", file=f)
            for gene in genes:
                print(gene.tag, gene.name, gene.location.start, gene.location.end, gene.strand, gene.score, gene.p_value, gene.protein_id, '"'+str(gene.product)+'"', sep=",", file=f)

    
    def refinment_it(self, p_value, k, N):
        """ Performs a p-value refinement iteration
            using a MM of order k, and a sampling size of N            
        """
        #instantiate MM generator
        generator = Generator(Generator.MCM)
        #number of pseudoreplicates with score > observed to reach p-value
        threshold = math.ceil(p_value * N)
        print(threshold)

        new_placements = []

        # for each window with a significant placement
        for placement in self.placements:
            # obtain the sequence (window) that generated the placement
            seq = placement.dna_sequence
            # generate N pseudoreplicates of the window sequence
            aux_seqs = (generator.generate(N, seq=seq, k=k))
            scores = []
            # perform placements on all pseudoreplicates
            for seq in aux_seqs:
                aux_placement = self.model.get_placement(seq)
                scores.append(aux_placement.energy)
            # sort scores and determine if placement score is above threshold
            scores.sort()
            score = scores[-threshold]
            if placement.energy >= score:
                new_placements.append(placement)

        #update list of placements with the ones making the cut
        self.placements = new_placements

    
    def refine(self, n=1, p_value=0.05, k=2):
        """ Performs iterative refinement of the p-value estimated for a placement
            using a sampling size of n (number of pseudosequences)
            a terminal p-value to be reached (p_value) and a Markov Model order (k).
        """
        # estimate the best MM order (k) to used based on user preference
        # and length of the sequence it is trained on [shorter sequences lead to
        # lower order Markov Models]
        aux_k = self.infer_k(self.placements[0].dna_sequence)
        k = min(k, aux_k)
        print("Refinement k", k)
        
        N = 10    # initial sampling size
        i=1       # iteration
        while N <= n:    # while not reaching terminal sampling size
            if config.VERBOSE:
                print("Refinement iteration", i)
            self.refinment_it(p_value, k, N)        #make refinement iteration
            if config.VERBOSE:
                print("Remaining placements:", len(self.placements))
            N*=10
            i+=1

    
    def get_threshold_score(self, WS, p_value=0.05, N=100):
        """ Computes initial score threshold for a given model.
            Estimates placement score on N random sequences and computes
            p-value from this, then determines
        """
        score = 0
        threshold = int(N*p_value)
        scores = []
        generator = Generator(Generator.RANDOM)
        seqs = generator.generate(N, n=WS)        # generate N random sequences

        # for each sequence, obtain optimal placement
        for seq in seqs:
            placement = self.model.get_placement(seq)
            scores.append(placement.energy)

        # sort scores and grab the score that corresponds to p-value limit
        scores.sort()
        self.threshold_scores = scores
        score = scores[-threshold]
        print("Threshold score", score)
        return score


    # def discard_placements(self, threshold):
    #     new_placements_forw = []
    #     new_placements_back = []

    #     for placement in self.placements_forw:
    #         if placement.energy >= threshold:
    #             new_placements_forw.append(placement)

    #     for placement in self.placements_back:
    #         if placement.energy >= threshold:
    #             new_placements_back.append(placement)
        
    #     self.placements_forw = new_placements_forw
    #     self.placements_back = new_placements_back


    def infer_k(self, seq):
        L = len(seq)
        return int(math.log(L/self.M, 4))-1
    
    def compute_p_value(self, score):
        count = len(list(filter(lambda x: x > score, self.threshold_scores)))

        return count/len(self.threshold_scores)
