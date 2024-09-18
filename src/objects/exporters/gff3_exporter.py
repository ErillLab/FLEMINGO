from .exporter import Exporter

class GFF3_Exporter(Exporter):

    def export(self, genome, placements: list) -> None:
        """ Exports the results in BED format. The fields dumped in
            this format are:

                1.- Sequence ID
                2.- Origen of the model
                3.- Type of element
                4.- Start position
                5.- End position
                6.- Score
                7.- Strand
                8.- Irrelevant
                9.- Additional attributes:
                    - Element ID
                    - Parent element ID
                    - Color
                    - Genes information
        """
        f = open(self.filename, "w")
        i = 0
        # dump placements to GFF3 file
        for placement in placements:
            # dump placement information as GFF3 element
            print(genome.name, "FLEMINGO", "TF_binding_site", placement.start_pos()+1, placement.end_pos(), placement.energy, placement.strand, 1, "ID=TFBS%d;height=7"%(i), sep="\t", file=f, end="")
            # dump left gene information
            l_gene = placement.l_gene
            if l_gene != None:
                print(";LG_NAME=%s"%l_gene.name, "LG_DISTANCE=%d"%(placement.end_pos() - l_gene.location.end), "LG_TAG=%s"%l_gene.tag, sep=";", file=f, end="")

            # dump right gene information
            r_gene = placement.r_gene
            if r_gene != None:
                print(";RG_NAME=%s"%r_gene.name, "RG_DISTANCE=%d"%(r_gene.location.start - placement.start_pos()), "RG_TAG=%s"%r_gene.tag, sep=";", file=f, end="")
            print("", file=f)

            # dump placement elements information as independent GFF3 element with the placement as parent
            for j in range(len(placement.recognizers_positions)):
                if j != 0:
                    connector = placement.connectors_positions[j-1]
                    print(genome.name, "FLEMINGO", "connector", connector[0]+1, connector[1], placement.connectors_scores[j-1], placement.strand, 1, "ID=TFBS%d_C%d;Parent=TFBS%d;height=1;color=rgba(0, 0, 0, 0)"%(i, j, i), sep="\t", file=f)
                recognizer = placement.recognizers_positions[j]
                # select color depending on element type
                if placement.recognizer_types[j] == "p":
                    elem_type = "PSSM"
                    color="rgb(55, 126, 184)"
                elif placement.recognizer_types[j] == "m":
                    elem_type = "mgw"
                    color="rgb(152, 78, 163)"
                elif placement.recognizer_types[j] == "t":
                    elem_type = "prot"
                    color="rgb(255, 127, 0)"
                elif placement.recognizer_types[j] == "h":
                    elem_type = "helt"
                    color="rgb(77, 175, 74)"
                elif placement.recognizer_types[j] == "r":
                    elem_type = "roll"
                    color="rgb(228, 26, 28)"
                print(genome.name, "FLEMINGO", elem_type, recognizer[0]+1, recognizer[1], placement.recognizers_scores[j], placement.strand, 1, "ID=TFBS%d_R%d;Parent=TFBS%d;height=7;color=%s"%(i, j, i, color), sep="\t", file=f)
            i+=1
        f.close()
