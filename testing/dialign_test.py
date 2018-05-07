import utilities


class PositionDict(object):
    """
    Contains information about a the positions and lengths of a sequences aligned columns

    """
    sequence_name = ""
    sequence_num = 0
    positions = []
    lengths = []

    def __init__(self, sequence_name, sequence_num, positions=None, lengths=None):
        self.sequence_name = sequence_name
        self.sequence_num = sequence_num
        self.positions = positions if type(positions) == list else []
        self.lengths = lengths if type(lengths) == list else []

    def add_entry(self, position, length):
        self.positions.append(position)
        self.lengths.append(length)


def create_dialign_anchor_file(alignment, outpath):
    """
    Take a structurally aligned alignment file and generate an anchor file of the aligned columns, to be used with
    DIALIGN to align the intermediate regions
    :param alignment: The structural alignment
    :param outpath: The path to write the anchor file to
    :return:
    """

    template = alignment[0].seq

    anchor_dict = get_anchor_dict(template)
    position_list = get_position_list(alignment, anchor_dict)

    seq_num = 0
    with open(outpath, "w+") as outfile:
        while seq_num + 1 < len(position_list):
            current_seq = position_list[seq_num]
            next_seq = position_list[seq_num + 1]
            for pos in range(0, len(current_seq.positions)):
                string = "%d %d %d %d %d %d \n" % (current_seq.sequence_num, next_seq.sequence_num,
                                                   current_seq.positions[pos], next_seq.positions[pos],
                                                   current_seq.lengths[pos], 100)
                outfile.write(string)
            seq_num += 2

        if len(position_list) % 2 != 0:
            current_seq = position_list[seq_num - 1]
            next_seq = position_list[seq_num]
            for pos in range(0, len(current_seq.positions)):
                string = "%d %d %d %d %d %d \n" % (current_seq.sequence_num, next_seq.sequence_num,
                                                   current_seq.positions[pos], next_seq.positions[pos],
                                                   current_seq.lengths[pos], 100)
                outfile.write(string)


def get_anchor_dict(template):
    """
    Take a singled aligned sequence and create a dictionary that maps the start position of a set of aligned columns to
    the length of the uninterrupted run of aligned columns
    :param template:
    :return:
    """
    anchor_dict = {}
    pos = 0
    while pos < len(template):
        length = 1
        if template[pos].isupper():
            start_pos = pos
            while pos + 1 < len(template) and template[pos + 1].isupper():
                pos += 1
                length += 1
            anchor_dict[start_pos] = length
        pos += 1
    return anchor_dict


def get_position_list(alignment, anchor_dict):
    """
    Take an alignment and a dictionary that maps the template anchors to their position with gaps and create a list of
    position dictionaries for each sequence that takes into account the positions once we remove gaps
    :param alignment:
    :param anchor_dict:
    :return:
    """
    position_list = []

    for num, sequence in enumerate(alignment):
        position_dict = PositionDict(sequence.name, num + 1)

        # Add in all the positions and lengths for the current sequence
        for pos, length in sorted(anchor_dict.items()):
            corrected_pos = pos - sequence.seq[0:pos].count("-") + 1
            position_dict.add_entry(corrected_pos, length)
        # Collate all the position dicts into the one list
        position_list.append(position_dict)

    return position_list


aln = utilities.load_alignment("dialign.aln")
create_dialign_anchor_file(aln, "./dialign.anc")
