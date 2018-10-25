import utilities

def correct_phylip_tree(phylip_file_path, phylip_dict_path, outpath):
    """
    Using a dictionary mapping full id to PHYLIP id, change a PHYLIP annotated tree back to the full ID
    :param records: The alignment we will retrieve the original IDs from
    :param phylip_file: The tree with the truncated PHYLIP ids
    :return:
    """
    # records = utilities.load_sequences(records_path)
    phylip_tree = utilities.load_tree(phylip_file_path)
    phylip_correction_dict = utilities.open_python_object(phylip_dict_path)


    for node in phylip_tree:
        if node.name in phylip_correction_dict:
            node.name = phylip_correction_dict[node.name]

    phylip_tree.write(outfile=outpath)


def generate_phylip_correction_dictionary(records, outpath=""):
    """
    Create a dictionary to map full id to a generated PHYLIP id, for cases when PHYLIP ids would be non-unique
    :return:
    """

    phylip_correction_dict = {}
    for record in records.values():
        if record.name[0:10] in phylip_correction_dict:
            phylip_correction_dict[utilities.random_string(10)] = record.name

        else:
            phylip_correction_dict[record.name[0:10]] = record.name

    if outpath:
        utilities.save_python_object(phylip_correction_dict, outpath)

    return phylip_correction_dict


def translate_to_phylip(records, phylip_dict):
    translated_records = {}
    for short, full in phylip_dict.items():
        translated_records[short] = records[full]
        translated_records[short].id = short
        translated_records[short].name = short
        translated_records[short].description = ""
    return translated_records


def write_out_phylip(records, outpath, outformat="fasta", dict_outpath = ""):
    phylip_correction_dict = generate_phylip_correction_dictionary(records, dict_outpath)
    translated_records = translate_to_phylip(records, phylip_correction_dict)
    out_records = map_dict_to_records(translated_records)
    SeqIO.write(out_records, outpath, outformat)