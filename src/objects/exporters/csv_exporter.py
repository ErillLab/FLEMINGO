from .exporter import Exporter

class CSV_Exporter(Exporter):

    def export(self, genome, placements: list) -> None:
        """ Exports the results in CSV format. The fields dumped in
            this format are:

                1.- Sequence ID
                2.- Placement ID
                3.- Start coordinate of the placement
                4.- End coordinate of the placement
                5.- Strand
                6.- Placement score
                7.- Placement p-value
                8..13.- Left gene information (Name, Tag, Distance, Strand, Protein Accession ID, Product)
                14..18- Right gene information (Name, Tag, Distance, Strand, Protein Accession ID, Product)
                19.- Placement sequence
                20..- Subsequences of placement elements
        """
        
        f = open(self.filename, "w")
        # dump header
        if len(placements) == 0:
            print("Sequence", "ID", "Start", "End", "Strand", "Score", "p-value", "LG_NAME", "LG_TAG", "LG_DISTANCE", "LG_STRAND", "LG_POTEIN_ID", "LG_PRDUCT", "RG_NAME", "RG_TAG", "RG_DISTANCE", "RG_STRAND", "RG_PROTEIN_ID", "LG_PRODUCT", "PLACEMENT_SEQUENCE", sep=",", file=f)
        else:
            print("Sequence", "ID", "Start", "End", "Strand", "Score", "p-value", "LG_NAME", "LG_TAG", "LG_DISTANCE", "LG_STRAND", "LG_POTEIN_ID", "LG_PRDUCT", "RG_NAME", "RG_TAG", "RG_DISTANCE", "RG_STRAND", "RG_PROTEIN_ID", "LG_PRODUCT", "PLACEMENT_SEQUENCE", self._get_model_elems(placements[0]), sep=",", file=f)
        i = 0

        # dump placements
        for placement in placements:
            # dump placement information
            print(genome.name, "TFBS"+str(i), placement.start_pos()+1, placement.end_pos(), placement.strand, placement.energy, placement.p_value, file=f, sep=",", end=",")
            # dump left gene information
            l_gene = placement.l_gene
            if l_gene != None:
                print(l_gene.name, l_gene.tag, (placement.end_pos() - l_gene.location.end), l_gene.strand, l_gene.protein_id, '"'+str(l_gene.product)+'"', sep=",", file=f, end=",")
            else: 
                print("-", "-", "-", "-", "-", "-", sep=",", file=f, end=",")

            # dump right gene information
            r_gene = placement.r_gene
            if r_gene != None:
                print(r_gene.name, r_gene.tag, (r_gene.location.start - placement.start_pos()), r_gene.strand, r_gene.protein_id, '"'+str(r_gene.product)+'"', sep=",", file=f, end=",")
            else: 
                print("-", "-", "-", "-", "-", "-", sep=",", file=f, end=",")

            # dump placement sequence and the subsequences for each model element
            print(",".join(self._report_placement(placement)), file=f)

            i+=1
        f.close()

    def _report_placement(self, placement) -> tuple:
        """ Returns two strings:
            - a concatenation of all the nucleotides bound by the chain structure
            - same but with a comma to separate between the subsequences bound by different recognizers
        """

        seq = placement.dna_sequence
        res = []
        offset = placement.window[0]

        
        for i in range(len(placement.connectors_positions)):
            # get recognizer subsequence
            ini = placement.recognizers_positions[i][0] - offset
            end = placement.recognizers_positions[i][1] - offset
            res.append(seq[ini:end])

            # get connector subsequence
            ini = placement.connectors_positions[i][0] - offset
            end = placement.connectors_positions[i][1] - offset
            res.append(seq[ini:end])

        # get last recognizer subsequence
        ini = placement.recognizers_positions[-1][0] - offset
        end = placement.end_pos() - offset
        res.append(seq[ini:end])
        return "".join(res), ",".join(res)
    
    def _get_model_elems(self, placement):
        """ Generates a coma separated list with the types of 
            model elements
        """
        res = ""
        
        for i in range(len(placement.recognizer_types)-1):
            if placement.recognizer_types[i] == "p":
                res += "PSSM"
            elif placement.recognizer_types[i] == "m":
                res += "MGW"
            elif placement.recognizer_types[i] == "t":
                res += "PORT"
            elif placement.recognizer_types[i] == "h":
                res += "HELT"
            elif placement.recognizer_types[i] == "r":
                res += "ROLL"
            res += ",CONNECTOR,"

        if placement.recognizer_types[-1] == "p":
            res += "PSSM"
        elif placement.recognizer_types[-1] == "m":
            res += "MGW"
        elif placement.recognizer_types[-1] == "t":
            res += "PORT"
        elif placement.recognizer_types[-1] == "h":
            res += "HELT"
        elif placement.recognizer_types[-1] == "r":
            res += "ROLL"
        return res
