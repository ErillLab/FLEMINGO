from .exporter import Exporter

class BED_Exporter(Exporter):
    def export(self, genome, placements: list) -> None:
        """ Exports the results in BED format. The fields dumped in
            this format are:

                1.- Sequence ID
                2.- Start coordinate of the placement
                3.- End coordinate of the placement
                4.- Placement ID
                5.- Placement score
                6.- Strand
                7.- Irrelevant
                8.- Irrelevant
                9.- Irrelevant
                10.- Number of recognizers
                11.- Recognizers lengths
                12.- Recognizer start position relative to placement start
        """
        f = open(self.filename, "w")
        i = 0
        for placement in placements:
            # dump the information of the placement
            print(genome.name, placement.start_pos(), placement.end_pos(), i, placement.energy, placement.strand, 0, 0, "0,0,255", len(placement.recognizers_positions), sep="\t", file=f, end="\t")
            # dump the size of the recognizers
            for element in placement.recognizers_positions[:-1]:
                size = element[1]-element[0]
                print(size, end=",", file=f)
            size = placement.end_pos()-placement.recognizers_positions[-1][0]
            print(size, end="\t", file=f)
            
            # dump the relative positions of the recognizers
            for element in placement.recognizers_positions[:-1]:
                print(int(element[0])-int(placement.start_pos()), end=",", file=f)
            print(int(placement.recognizers_positions[-1][0])-int(placement.start_pos()), end="\n", file=f)
            i+=1
        f.close()
