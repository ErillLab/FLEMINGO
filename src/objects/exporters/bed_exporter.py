from .exporter import Exporter

class BED_Exporter(Exporter):

    def export(self, genome, placements: list) -> None:
        f = open(self.filename, "w")
        i = 0
        for placement in placements:
            print(genome.name, placement.recognizers_positions[0][0], placement.recognizers_positions[-1][1], i, placement.energy, placement.strand, 0, 0, "0,0,255", len(placement.recognizers_positions), sep="\t", file=f, end="\t")
            for element in placement.recognizers_positions[:-1]:
                size = element[1]-element[0]
                print(size, end=",", file=f)
            size = placement.recognizers_positions[-1][1]-placement.recognizers_positions[-1][0]
            print(size, end="\t", file=f)
            
            for element in placement.recognizers_positions[:-1]:
                print(int(element[0])-int(placement.recognizers_positions[0][0]), end=",", file=f)
            print(int(placement.recognizers_positions[-1][0])-int(placement.recognizers_positions[0][0]), end="\n", file=f)
            i+=1
        f.close()
