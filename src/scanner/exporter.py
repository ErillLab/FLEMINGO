import FMS.utils as utils
import json
import FMS.config as config

global results
results = {}
placements = []

def export_results(GL, model, WS, inter_reg,):
    results["genome_length"] = GL
    results["window_size"] = WS
    results["CI"] = config.LENGTH_CI
    results["format"] = config.INPUT_TYPE
    results["model"] = model.to_json()
    # results["organism"] = orgnfact.export_organisms([organism])
    # if config.SCAN_ONLY_INTERGENIC:
    #     results["intergenic_regions"] = inter_reg
    results["placements"] = placements
    with open(config.OUTPUT_FILE, "w") as f:
        json.dump(results, f, indent=4)

def export_placement(wstart, wend, placement, ln_gene, rn_gene, strand="+"):
    placement_json = {}
    placement_json["window_start"] = wstart
    placement_json["window_end"] = wend
    placement_json["strand"] = strand

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

    placements.append(placement_json)
