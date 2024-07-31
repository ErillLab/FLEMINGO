import FMS.utils as utils
from FMS.objects.genbank.genbank_object import GenBank
import FMS.config as config

class Genome:

    def __init__(self, filename, format="fasta"):
        self.genes = []
        self.inter_regions = []
        self.upstream_regions = []
        # self.cds_forw = []
        # self.cds_back = []
        self.cds = []
        self.filename = filename

        if format == "fasta":
            self._fasta_genome()
        elif format == "genbank":
            self._genbank_genome()

        self.length = len(self.seq)
        if config.SCAN_MODE == "INTERGENIC":
            self.inter_regions = self._calc_intergenic_regions(self.cds, self.length)
        elif config.SCAN_MODE == "UPSTREAM_REGIONS":
            self.upstream_regions = self._calc_upstream_regions(self.genes, config.UPSTREAM_REGION_SIZE)
        # self.inter_regions_forw = self._calc_intergenic_regions(self.cds_forw, self.length)
        # self.inter_regions_back = self._calc_intergenic_regions(self.cds_back, self.length)

    def _fasta_genome(self):
        self.seq, self.name = utils.read_fasta_file(self.filename)[0]


    def _genbank_genome(self):
        genbank = GenBank(self.filename)
        self.seq = genbank.seqs[0].seq
        self.genes = genbank.seqs[0].genes
        self.cds = genbank.seqs[0].cds
        self.cds_forw = genbank.seqs[0].cds_forw
        self.cds_back = genbank.seqs[0].cds_back
        self.name = (genbank.seqs[0].name)

    def _calc_intergenic_regions(self, cds_a, N):
        res = [[0, N]]
        for cds in cds_a:
            ini = int(cds.location.start)
            end = int(cds.location.end)
            for pos in res:
                if ini > pos[0] and ini < pos[1]:
                    if end < pos[1]:
                        res.append([end, pos[1]])
                        pos[1] = ini
                        break
                    pos[1] = ini
                if end < pos[1] and end > pos[0]:
                    if ini > pos[0]:
                        res.append([pos[0], ini])
                    pos[0] = end
                    break
        return res
    
    def _calc_upstream_regions(self, genes, N):
        ups_regs = []
        for gene in genes:
            ini = int(gene.location.start)
            end = int(gene.location.end)
            if gene.strand == "+":
                ups_regs.append([end, min(end+N, self.length)])
            else:
                ups_regs.append([max(0, ini-N), ini])
        ups_regs.sort(key=lambda x: x[0])
    
        res = []
        aux = ups_regs[0]
    
        for reg in ups_regs[1:]:
            if reg[0] <= aux[1]:
                aux = (aux[0], max(aux[1], reg[1]))
            else:
                res.append(aux)
                aux = reg
        res.append(aux)
        return res

    def get_nearest_genes(self, start, end, strand="+"):
        ln_gene_diff = self.length
        rn_gene_diff = self.length
        ln_gene = None
        rn_gene = None
        for gene in self.genes:
            # if gene.strand != strand:
            #     continue
            if gene.location.start - start < rn_gene_diff and gene.location.start - start > 0:
                rn_gene_diff = gene.location.start - start
                rn_gene = gene
            if end - gene.location.end < ln_gene_diff and end - gene.location.end > 0:
                ln_gene_diff = end - gene.location.end
                ln_gene = gene

        if rn_gene != None and rn_gene.strand != "+":
            rn_gene = None
        if ln_gene != None and ln_gene.strand != "-":
            ln_gene = None
        return ln_gene, rn_gene
