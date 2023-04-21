import objects.organism_factory as organism_factory
import json
import models.models as models  

   
def test_placement(organism, sequence):
    placement = organism.get_placement(sequence)
    placement.print_placement(stdout=True)
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
    org_list = org_fac.import_organisms("organism.json")
    return org_list[0]


def main():
    org = build_organism()
    org.flatten()
    placement = test_placement(org, "AAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTTAAGGTT")
    try:
        placement = test_placement(org, "AAGG")
    except:
        pass
    return

if __name__ == '__main__':
    main()