from .exporter import Exporter

class CSV_Exporter(Exporter):

    def export(self, genome, placements: list) -> None:
        
        f = open(self.filename, "w")
        if len(placements) == 0:
            
            print("Sequence", "ID", "Start", "End", "Strand", "Score", "p-value", "LG_NAME", "LG_TAG", "LG_DISTANCE", "LG_STRAND", "LG_POTEIN_ID", "LG_PRDUCT", "RG_NAME", "RG_TAG", "RG_DISTANCE", "RG_STRAND", "RG_PROTEIN_ID", "LG_PRODUCT", "PLACEMENT_SEQUENCE", sep=",", file=f)
        else:
            print("Sequence", "ID", "Start", "End", "Strand", "Score", "p-value", "LG_NAME", "LG_TAG", "LG_DISTANCE", "LG_STRAND", "LG_POTEIN_ID", "LG_PRDUCT", "RG_NAME", "RG_TAG", "RG_DISTANCE", "RG_STRAND", "RG_PROTEIN_ID", "LG_PRODUCT", "PLACEMENT_SEQUENCE", self._get_model_elems(placements[0]), sep=",", file=f)
        i = 0
        for placement in placements:
            print(genome.name, "TFBS"+str(i), placement.recognizers_positions[0][0]+1, placement.recognizers_positions[-1][1], placement.strand, placement.energy, placement.p_value, file=f, sep=",", end=",")
            l_gene = placement.l_gene
            if l_gene != None:
                print(l_gene.name, l_gene.tag, (placement.recognizers_positions[-1][1] - l_gene.location.end), l_gene.strand, l_gene.protein_id, '"'+str(l_gene.product)+'"', sep=",", file=f, end=",")
            else: 
                print("-", "-", "-", "-", "-", "-", sep=",", file=f, end=",")

            r_gene = placement.r_gene
            if r_gene != None:
                print(r_gene.name, r_gene.tag, (r_gene.location.start - placement.recognizers_positions[0][0]), r_gene.strand, r_gene.protein_id, '"'+str(r_gene.product)+'"', sep=",", file=f, end=",")
            else: 
                print("-", "-", "-", "-", "-", "-", sep=",", file=f, end=",")

            print(",".join(self._report_placement(placement)), file=f)

            i+=1
        f.close()

    def _report_placement(self, placement):

        seq = placement.dna_sequence
        res = []
        offset = placement.window[0]

        for i in range(len(placement.connectors_positions)):
            ini = placement.recognizers_positions[i][0] - offset
            end = placement.recognizers_positions[i][1] - offset
            res.append(seq[ini:end])

            ini = placement.connectors_positions[i][0] - offset
            end = placement.connectors_positions[i][1] - offset
            res.append(seq[ini:end])

        ini = placement.recognizers_positions[-1][0] - offset
        end = placement.recognizers_positions[-1][1] - offset
        res.append(seq[ini:end])
        return "".join(res), ",".join(res)
    
    def _get_model_elems(self, placement):
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
