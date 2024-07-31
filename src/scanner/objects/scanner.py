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
        config.set_up(cfg_file=cfg_file)

    def load(self):
        self.genome = Genome(config.INPUT_GENOME_FILE, format=config.INPUT_TYPE)

        self.model = ChainModel.import_model(config.INPUT_ORGANISM_FILE)
        if config.VERBOSE:
            self.model.print()

        self.WS = self.model.get_max_length(config.LENGTH_CI)
        print(self.WS)

        self.threshold_p_value = config.THRESHOLD_P_VALUE
        self.sample_size = config.SAMPLE_SIZE
        self.threshold = self.get_threshold_score(self.WS, N=self.sample_size, p_value=self.threshold_p_value)
 
    def scan(self):
        if config.SCAN_MODE == "INTERGENIC":
            regions = self.genome.inter_regions
        
        elif config.SCAN_MODE == "UPSTREAM_REGIONS":
            regions = self.genome.upstream_regions
        else:
            regions = [[0, self.genome.length]]

        # prev_pos_score = 0
        # prev_pos_loc = 0
        # prev_neg_score = 0
        # prev_neg_loc = 0
        prev_score = -1
        prev_loc = -1

        for reg in regions:

            size = reg[1]-reg[0]
            aux_N = ((size // self.WS)+1)*self.WS
            offset = reg[0] + size//2 - aux_N//2

            if offset < 0:
                offset = 0
            if aux_N > self.genome.length:
                aux_N = self.genome.length
            for n in range(0, aux_N - (self.WS >> 1), self.WS>>1):
                window = str(self.genome.seq[n+offset:n+offset+self.WS])
                if len(window) < self.WS>>1:
                    break
                for_placement = self.model.get_placement(window)
                if config.SCAN_REVERSE:
                    window = str(Seq(window).reverse_complement())
                    rev_placement = self.model.get_placement(window)

                if not config.SCAN_REVERSE or for_placement.energy >= rev_placement.energy:
                    placement = for_placement
                    placement.strand = "+"
                else:
                    placement = rev_placement
                    placement.strand = "-"
                    # prev_neg_loc = offset + n + self.WS - rev_placement.recognizers_positions[0][0] 
                    # prev_neg_score = rev_placement.energy
                    placement = utils.reverse_placement(placement, len(window))

                if placement.energy < self.threshold or (prev_score == placement.energy and prev_loc == offset + n + for_placement.recognizers_positions[0][0]):
                    continue
                        

                
                # if for_placement.energy >= self.threshold and (prev_pos_loc != offset + n + for_placement.recognizers_positions[0][0] or prev_pos_score != for_placement.energy):
                #     prev_pos_loc = offset + n + for_placement.recognizers_positions[0][0]
                #     prev_pos_score = for_placement.energy
                    # l_gene, r_gene = self.genome.get_nearest_genes(offset + n + for_placement.recognizers_positions[0][0], offset + n + for_placement.recognizers_positions[-1][1])
                    # for_placement.r_gene = r_gene
                    # for_placement.l_gene = l_gene
                    # if config.VERBOSE:
                    #     utils.print_placement_info(for_placement, l_gene, r_gene, n + offset)
                    # for_placement.add_offset(n + offset)
                    # for_placement.window=[n+offset, n+offset+self.WS]
                    # for_placement.p_value = self.compute_p_value(for_placement.energy)
                    # for_placement.strand = "+"
                    # self.placements_forw.append(for_placement)

                # if config.SCAN_REVERSE:
                    # window = str(Seq(window).reverse_complement())
                    # rev_placement = self.model.get_placement(window)
                    # if rev_placement.energy >= self.threshold and (prev_neg_loc != offset + n + self.WS - rev_placement.recognizers_positions[0][0] or prev_neg_score != rev_placement.energy):
                    #     prev_neg_loc = offset + n + self.WS - rev_placement.recognizers_positions[0][0] 
                    #     prev_neg_score = rev_placement.energy
                        # rev_placement = utils.reverse_placement(rev_placement, len(window))
                        # l_gene, r_gene = self.genome.get_nearest_genes(offset + n + rev_placement.recognizers_positions[0][0], offset + n + rev_placement.recognizers_positions[-1][1], strand="-")
                        # rev_placement.r_gene = r_gene
                        # rev_placement.l_gene = l_gene
                        # if config.VERBOSE:
                        #     utils.print_placement_info(rev_placement, l_gene, r_gene, n + offset)
                        # rev_placement.add_offset(n + offset)
                        # rev_placement.window=[n+offset, n+offset+self.WS]
                        # rev_placement.p_value = self.compute_p_value(rev_placement.energy)
                        # rev_placement.strand = "-"
                        # self.placements_back.append(rev_placement)
                l_gene, r_gene = self.genome.get_nearest_genes(offset + n + placement.recognizers_positions[0][0], offset + n + placement.recognizers_positions[-1][1], strand=placement.strand)
                placement.r_gene = r_gene
                placement.l_gene = l_gene
                if config.VERBOSE:
                    utils.print_placement_info(placement, l_gene, r_gene, n + offset)
                placement.add_offset(n + offset)
                placement.window=[n+offset, n+offset+self.WS]
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
        generator = Generator(Generator.MCM)
        threshold = math.ceil(p_value * N)
        print(threshold)

        new_placements = []
        # new_placements_forw = []
        # new_placements_back = []

        for placement in self.placements:
            seq = placement.dna_sequence
            aux_seqs = (generator.generate(N, seq=seq, k=k))
            scores = []
            for seq in aux_seqs:
                aux_placement = self.model.get_placement(seq)
                scores.append(aux_placement.energy)
            scores.sort()
            score = scores[-threshold]
            if placement.energy >= score:
                new_placements.append(placement)


        # for placement in self.placements_back:
        #     seq = placement.dna_sequence[::-1]
        #     aux_seqs = (generator.generate(N, seq=seq, k=k))
        #     scores = []
        #     for seq in aux_seqs:
        #         aux_placement = self.model.get_placement(seq)
        #         scores.append(aux_placement.energy)
        #     scores.sort()
        #     score = scores[-threshold]
        #     if placement.energy >= score:
        #         new_placements_back.append(placement)
        self.placements = new_placements
        # self.placements_forw = new_placements_forw
        # self.placements_back = new_placements_back

    def refine(self, n=1, p_value=0.05, k=2):
        aux_k = self.infer_k(self.placements[0].dna_sequence)
        # aux_k = self.infer_k(self.placements_forw[0].dna_sequence)
        k = min(k, aux_k)
        print("Refinment k", k)
        N = 10
        i=1
        while N <= n:
            if config.VERBOSE:
                print("Refinment iteration", i)
                self.refinment_it(p_value, k, N)
            if config.VERBOSE:
                print("Remaining placements:", len(self.placements))
                # print("Remaining placements:", len(self.placements_forw) + len(self.placements_back))
            N*=10
            i+=1

    def get_threshold_score(self, WS, p_value=0.05, N=100):
        score = 0
        threshold = int(N*p_value)
        scores = []
        generator = Generator(Generator.RANDOM)
        seqs = generator.generate(N, n=WS)

        for seq in seqs:
            placement = self.model.get_placement(seq)
            scores.append(placement.energy)
        
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