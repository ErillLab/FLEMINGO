from .exporter import Exporter
import json
from objects import config

class JSON_Exporter(Exporter):

    def export(self, genome, placements: list, model) -> None:
        self.results = {}
        self.placements = []

        for placement in placements:
            self.placement_to_json(placement)

        self.results["genome_length"] = genome.length
        self.results["CI"] = config.LENGTH_CI
        self.results["format"] = config.INPUT_TYPE
        self.results["model"] = model.to_json()
        self.results["placements"] = self.placements
        with open(self.filename, "w") as f:
            json.dump(self.results, f, indent=4)

    def placement_to_json(self, placement: list) -> dict:
        ln_gene = placement.l_gene
        rn_gene = placement.r_gene

        placement_json = {}
        placement_json["window_start"] = placement.window[0]
        placement_json["window_end"] = placement.window[1]
        placement_json["strand"] = placement.strand

        placement_json["placement_start"] = placement.recognizers_positions[0][0]
        placement_json["placement_end"] = placement.recognizers_positions[-1][1]
        placement_json["length"] = placement.recognizers_positions[-1][1] - placement.recognizers_positions[0][0]
        placement_json["score"] = placement.energy

        if config.INPUT_TYPE == "genbank":
            left_gene = None
            right_gene = None

            if ln_gene != None:
                left_gene = {}
                left_gene["name"] = ln_gene.name
                left_gene["distance"] = placement.recognizers_positions[0][0] - ln_gene.location.end
                left_gene["strand"] = "+" if (ln_gene.location.strand==1) else "-"
                left_gene["gene_synonym"] = ln_gene.synonym
                left_gene["tag"] = ln_gene.tag

            if rn_gene != None:
                right_gene = {}
                right_gene["name"] = rn_gene.name
                right_gene["distance"] = rn_gene.location.start - placement.recognizers_positions[-1][1]
                right_gene["strand"] = "+" if (rn_gene.location.strand==1) else "-"
                right_gene["gene_synonym"] = rn_gene.synonym
                right_gene["tag"] = rn_gene.tag

            placement_json["left_gene"] = left_gene
            placement_json["right_gene"] = right_gene

        placement_json["organism"] = []
        recognizer = {}
        recognizer["type"] = "Recognizer"
        recognizer["start"] = placement.recognizers_positions[0][0]
        recognizer["end"] = placement.recognizers_positions[0][1]
        recognizer["size"] = placement.recognizers_positions[0][1] - placement.recognizers_positions[0][0]
        recognizer["energy"] = placement.recognizers_scores[0]
        placement_json["organism"].append(recognizer)

        for i in range(len(placement.connectors_positions)):
            connector = {}
            connector["type"] = "Connector"
            connector["start"] = placement.connectors_positions[i][0]
            connector["end"] = placement.connectors_positions[i][1]
            connector["size"] = placement.connectors_positions[i][1] - placement.connectors_positions[i][0]
            connector["energy"] = placement.connectors_scores[i]
            placement_json["organism"].append(connector)

            recognizer = {}
            recognizer["type"] = "Recognizer"
            recognizer["start"] = placement.recognizers_positions[i+1][0]
            recognizer["end"] = placement.recognizers_positions[i+1][1]
            recognizer["size"] = placement.recognizers_positions[i+1][1] - placement.recognizers_positions[i+1][0]
            recognizer["energy"] = placement.recognizers_scores[i+1]
            placement_json["organism"].append(recognizer)

        self.placements.append(placement_json)
