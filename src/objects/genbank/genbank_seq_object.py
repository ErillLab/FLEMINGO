from .genbank_gene_object import GenBankGene
from .genbank_cds_object import GenBankCDS

class GenBankSeq:

    def __init__(self, info):
        """ Initializes GenBank sequence object
        """
        self.genes = []
        self.cds_back = []
        self.cds_forw = []
        self.cds = []
        self._parse_genbank_info(info)

    def _parse_genbank_info(self, info):
        """ Gets some attributes from the information 
        getted from the genbank file
        """
        self._parse_features(info.features)
        self.seq = str(info.seq)
        self.name = info.id

    def _parse_features(self, features):
        """ Parses the "features" of the genbank file,
            specifically genes and coding regions
        """
        for feature in features:
            # creaste a GenBank gene object in case of being a gene feature
            if feature.type == "gene":
                gene = GenBankGene(feature)
                self.genes.append(gene)
            # creaste a GenBank cds object in case of being a cds feature
            elif feature.type == "CDS":
                cds = GenBankCDS(feature)
                self.cds.append(cds)
                if cds.tag == gene.tag:
                    gene.product = cds.product
                    gene.protein_id = cds.protein_id
