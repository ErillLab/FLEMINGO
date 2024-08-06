from .exporter import Exporter
import json
from objects import config

class JSON_Exporter(Exporter):

    def export(self, genome, placements: list, model) -> None:
        """ Exports the results in JSON format. The fields dumped in
            this format are:

                1.- Genome length
                2.- Confidence interval used to determine the window size
                3.- Input genome format
                4.- Chain model used for placement detection
                5.- Significative placements
        """
        self.results = {}
        self.placements = []

        # generate a list of dictionaries with placements information
        for placement in placements:
            self.placement_to_json(placement)

        self.results["genome_length"] = genome.length
        self.results["CI"] = config.LENGTH_CI
        self.results["format"] = config.INPUT_TYPE
        self.results["model"] = model.to_json()
        self.results["placements"] = self.placements

        # dump generated dictionary to json file
        with open(self.filename, "w") as f:
            json.dump(self.results, f, indent=4)

    def placement_to_json(self, placement: list) -> dict:
        """ Creates a dictionary with the placement information
            to later dump it to json file

            placement: {
                window start,
                window end,
                strand,
                placement start,
                placement end,
                length,
                score,
                left gene information, (only GENBANK)
                right gene information, (only GENBANK)
                chain model elements information
            }

            gene: {
                name,
                tag,
                distance,
                strand,
                gene synonym
            }

            element: {
                element type,
                element start,
                element end,
                score,
                length
            }
        """
        ln_gene = placement.l_gene
        rn_gene = placement.r_gene

        # add placement information to the dictionary
        placement_json = {}
        placement_json["window_start"] = placement.window[0]
        placement_json["window_end"] = placement.window[1]
        placement_json["strand"] = placement.strand

        placement_json["placement_start"] = placement.start_pos()
        placement_json["placement_end"] = placement.end_pos()
        placement_json["length"] = placement.end_pos() - placement.start_pos()
        placement_json["score"] = placement.energy

        # in case of having a genome in genbank format dump genes information
        if config.INPUT_TYPE == "genbank":
            left_gene = None
            right_gene = None

            if ln_gene != None:
                left_gene = {}
                left_gene["name"] = ln_gene.name
                left_gene["distance"] = placement.start_pos() - ln_gene.location.end
                left_gene["strand"] = "+" if (ln_gene.location.strand==1) else "-"
                left_gene["gene_synonym"] = ln_gene.synonym
                left_gene["tag"] = ln_gene.tag

            if rn_gene != None:
                right_gene = {}
                right_gene["name"] = rn_gene.name
                right_gene["distance"] = rn_gene.location.start - placement.end_pos()
                right_gene["strand"] = "+" if (rn_gene.location.strand==1) else "-"
                right_gene["gene_synonym"] = rn_gene.synonym
                right_gene["tag"] = rn_gene.tag

            placement_json["left_gene"] = left_gene
            placement_json["right_gene"] = right_gene

        # add to placement information dictionary, subdictionaries with placement elements information
        placement_json["chain"] = []
        recognizer = {}
        recognizer["type"] = "Recognizer"
        recognizer["start"] = placement.start_pos()
        recognizer["end"] = placement.recognizers_positions[0][1]
        recognizer["size"] = placement.recognizers_positions[0][1] - placement.start_pos()
        recognizer["energy"] = placement.recognizers_scores[0]
        placement_json["chain"].append(recognizer)

        for i in range(len(placement.connectors_positions)):
            connector = {}
            connector["type"] = "Connector"
            connector["start"] = placement.connectors_positions[i][0]
            connector["end"] = placement.connectors_positions[i][1]
            connector["size"] = placement.connectors_positions[i][1] - placement.connectors_positions[i][0]
            connector["energy"] = placement.connectors_scores[i]
            placement_json["chain"].append(connector)

            recognizer = {}
            recognizer["type"] = "Recognizer"
            recognizer["start"] = placement.recognizers_positions[i+1][0]
            recognizer["end"] = placement.recognizers_positions[i+1][1]
            recognizer["size"] = placement.recognizers_positions[i+1][1] - placement.recognizers_positions[i+1][0]
            recognizer["energy"] = placement.recognizers_scores[i+1]
            placement_json["chain"].append(recognizer)

        self.placements.append(placement_json)
