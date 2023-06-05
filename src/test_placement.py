import objects.organism_factory as organism_factory
import json
import models.models as models  

   
def test_placement(organism, sequence):
    print("")
    print("---------------------TESTING---------------------")
    print("Num Recognizers:", len(organism.recognizers))
    print("Sequence Length:", len(sequence))
    print("-----------------    RESULTS    -----------------")
    placement = organism.get_placement(sequence)
    print("-------------------------------------------------")
    print("")
    #placement.print_placement(stdout=True)
    return placement

def build_organism():
    config = []
    with open("config.json") as json_content:
        config = json.load(json_content)

    configOrganism = config["organism"]
    configOrganismFactory = config["organismFactory"]
    configConnector = config["connector"]
    configPssm = config["pssm"]
    configShape = config["shape"]
    
    org_fac = organism_factory.OrganismFactory( configOrganism, configOrganismFactory, configConnector, configPssm, 0, configShape )
    org_list = org_fac.import_organisms("organism2.json")
    return org_list[0]


def main():
    org = build_organism()
    print("ONLY LOOK AT AFTER THIS STATEMENT")
    org.flatten()
    scores1 = org.connectors_scores_flat
    placement = test_placement(org, "gtacaacatg" * 3)
    placement.print_placement(stdout=True)
    exit()
    org.set_pf_auc()
    scores2 = org.connectors_scores_flat
    placement = test_placement(org, "gtacaacatg" * 10)
    placement.print_placement(stdout=True)
    #placement = test_placement(org, "gtacaacatg" * 1000)

    return

if __name__ == '__main__':
    main()