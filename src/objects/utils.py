from Bio import SeqIO
import json
import os
import objects.config as config

def read_fasta_file(filename: str) -> list:
    """Reads a fasta file and returns an array of DNA sequences (strings)

    TODO: probably it can be useful to create our own Sequence object that
    creates the string and stores some properties from fasta format. Also
    we can adapt the current program to use Biopythons's Seq object.

    Args:
    filename: Name of the file that contains FASTA format sequences to read

    Returns:
    The set of sequences in string format
    """""
    dataset = []
    fasta_sequences = SeqIO.parse(open(filename), "fasta")
    for fasta in fasta_sequences:
        dataset.append([str(fasta.seq), str(fasta.name)])
        return dataset

def read_genbank_file(filename: str) -> list:
    """ Reads a GenBank file and returns an array of DNA sequences (strings)

    Args:
        filename: Name of the file that contains GenBank format sequences to read

    Returns:
        The set of sequences in string format
    """
    dataset = []
    genbank_sequences = SeqIO.parse(open(filename, "r"), "genbank")
    for genbank in genbank_sequences:
        dataset.append(genbank.seq)
        return dataset


def read_json_file(filename: str) -> dict:
    """Reads a JSON file and returns a dictionary with the content.

    Args:
        filename: Name of the json file to read
    Returns:
        Dictionary with the json file info
    """""
    with open(filename) as json_content:
        return json.load(json_content)

def i_am_main_process() -> bool:
    ''' Returns True if rank is None (which happens when the run is serial) or
    when rank is 0 (which happens when the run is parallel and the program is
    executed by the process with rank 0). '''
    return not config.rank

def check_dir(dir_path: str) -> None:
    ''' If the directory at the specified path doesn't exists, it's created.
    Any missing parent directory is also created. If the directory already
    exists, it is left un modified. This function works even when executed in
    parallel my several processes (makedir would crash if they try to make the
    directory at about the same time).
    '''
    os.makedirs(dir_path, exist_ok=True)

def print_ln(string, filepath=None, to_stdout=True) -> None:
    """Shows the string on stdout and write it to a file
    (like the python's logging modules does)

    Args:
        string: Information to print on stdout and file
        filepath: path to the file to export the string
    """

    if to_stdout:
        print(string)

    if filepath != None:
        # Here we are sure file exists
        _f = open(filepath, "a+")
        _f.write(string + "\n")
        _f.close()

def get_file_type(filename: str) -> str:
    """ Gets the file format from the filename extension
    """
    ext = filename.split(".")[-1]
    if ext == "fna" or ext == "fasta" or ext == "fas":
        return "fasta"
    elif ext == "gbff" or ext == "gb":
        return "genbank"

def reverse_placement(placement, WS: int):
    """ Reverse the placement coordinates. It is used to update 
        coordinates for placement on the reverse strand, switching
        from using reverse strand coordinates to forward strand 
        coordinates 
    """

    # reverse the recognizers coordinates and scores
    for position in placement.recognizers_positions:
        position[0], position[1] = WS - position[1], WS - position[0]
    placement.recognizers_positions = placement.recognizers_positions[::-1]
    placement.recognizers_scores = placement.recognizers_scores[::-1]

    # reverse the connectors coordinates and scores
    for position in placement.connectors_positions:
        position[0], position[1] = WS - position[1], WS - position[0]
    placement.connectors_positions = placement.connectors_positions[::-1]
    placement.connectors_scores = placement.connectors_scores[::-1]

    # reverse the connector lis type
    placement.recognizer_types = placement.recognizer_types[::-1]

    # reverse the window sequence
    placement.dna_sequence = placement.dna_sequence[::-1]
    return placement

def print_placement_info(placement) -> None:
    """ Prints the placement, and the nearest genes in
        case of having a GenBank file
    """

    placement.print_placement(stdout=True)
    if config.INPUT_TYPE == "genbank":
        if placement.l_gene != None:
            print("Left gene (%d):"%(placement.start_pos() + placement.window[0] - placement.l_gene.location.end))
            placement.l_gene.print()
        if placement.r_gene != None:
            print("Right gene (%d):"%(placement.r_gene.location.start - placement.end_pos() - placement.window[0]))
            placement.r_gene.print()

            