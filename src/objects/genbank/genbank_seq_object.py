from Bio import SeqIO
from .genbank_gene_object import GenBankGene
from .genbank_cds_object import GenBankCDS

class GenBankSeq:

    def __init__(self, info):
        self.genes = []
        self.cds_back = []
        self.cds_forw = []
        self.cds = []
        self._parse_genbank_info(info)

    def _parse_genbank_info(self, info):
        self._parse_features(info.features)
        self.seq = str(info.seq)
        self.name = info.id

    def _parse_features(self, features):
        for feature in features:
            if feature.type == "gene":
                gene = GenBankGene(feature)
                self.genes.append(gene)
            elif feature.type == "CDS":
                cds = GenBankCDS(feature)
                # if cds.strand == "+":
                #     self.cds_forw.append(cds)
                # else:
                #     self.cds_back.append(cds)
                self.cds.append(cds)
                # print(cds.tag, gene.tag)
                if cds.tag == gene.tag:
                    gene.product = cds.product
                    gene.protein_id = cds.protein_id

        # exit()
