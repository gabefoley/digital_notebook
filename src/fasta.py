from Bio import SeqIO, Entrez
from collections import defaultdict
import plotly
import utilities
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.graph_objs import *
plotly.tools.set_credentials_file(username='gabefoley', api_key='xS8qT0kIbIKDWt0BalOd')


def print_record_overview(records):
    """
    Print the total number of sequences and the average length of sequences
    :param records: The records to print information about
    """
    total_len = 0
    for record in records.values():
        total_len += len(record)
    print("There are %s sequences and the average length of sequence is %d" % (len(records), total_len / len(records)))


def build_species_count(*exclude_on, records, length=0):
    """
    Builds up a dictionary with species as keys and lists of record ids that are from that species as values
    :param exclude_on:
    :param records:
    :param length:
    :return:
    """
    species_count = defaultdict(list)
    for record in records.values():
        if len(str(record.seq)) > length and not any(exclude in record.description for exclude in exclude_on):
            if "[" in record.description:
                name = record.description.split(">")[0].split("[")[1].split("]")[0]
            else:
                name = "_".join(record.description.split("|")[1].split("_")[0:2])
            species_count[name].append(record.id)
    return species_count


def subset_records(*header_terms, records, length=0, mode="exclude"):
    """
    Return a subset of records, by defining terms that should or shouldn't be in the header and a minimum length

    :param header_terms: The terms to look for in the header
    :param records: The records to make a subset from
    :param length: The minimum length of sequences to keep
    :param mode: Whether to exclude or include on the presence of terms in the header
    :return:
    """
    new_records = {}
    for record in records.values():
        if mode == "exclude":
            if len(str(record.seq)) > length and not any(term in record.description for term in header_terms):
                new_records[record.name] = record
        elif mode == "include":
            if len(str(record.seq)) > length and any(term in record.description for term in header_terms):
                new_records[record.name] = record

    return new_records


def add_species_from_header(records):
    """ Add species names that exist within the header of the seq ID i.e. like [Amazona aestiva] to a dictionary of
    seqRecords
    :param records: The dictionary to update with species names
    :return:
    """

    for record in records:
        if "[" in records[record].description:
            species_name = records[record].description.split(">")[0].split("[")[1].split("]")[0]
            records[record].annotations["Species"] = species_name
    return records


def add_species_from_last_position(records):
    """ Add species names that exist within the last positiong of the seq ID to a dictionary of seqRecords
    :param records: The dictionary to update with species names
    :return:
    """

    for record in records:
        if "_" in records[record].description:
            species_name = records[record].description.split("_")[-1]
            records[record].annotations["Species"] = species_name
    return records


def get_records_without_species(records):
    """
    Get the entries in a seq records file that don't have a species annotated with them
    :param records: The records to check
    :return: The list of seq IDs missing species information
    """
    missing = []

    for record in records.values():
        if "Species" not in record.annotations:
            missing.append(record.id)

    return missing


def report_records_without_species(records):
    """
    Print out statements about whether entries in a seq record file have annotated species
    :param records:
    :return:
    """
    missing = get_records_without_species(records)
    if missing:
        print("The following records do not have species annotated", [x for x in missing])
    else:
        print("All of the records have a species annotated")


def check_cd_hit_output(records, cdhit_output):
    """
    Reports on the species information that occurs within clusters in CD-hit output
    :param records:
    :param cdhit_output:
    :return:
    """

    if get_records_without_species(records):
        records = add_species_from_header(records)
        if get_records_without_species(records):
            print("Some of your entries are missing their species annotation")
            return

    clusters = defaultdict(list)
    current = None

    # Read in and organise the CD-HIT output

    for line in cdhit_output:
        if "Cluster" in line:
            current = line

        if "Cluster" not in line:

            clusters[current].append(re.sub(r"\s+", " ", line))

    # Only interested in clusters with more than one sequence
    for cluster in clusters.values():
        if len(cluster) > 1:
            seq_list = []
            for seq in cluster:
                trim = re.search('>(.*)\.\.\.', seq)
                seq_id = trim.group(1).strip()
                if "*" in seq:
                    seq_info = "Retained " + seq_id + " " + records[seq_id].annotations["Species"]
                    seq_list.insert(0, seq_info)  # Add the retained sequence to the start of the seq list
                else:
                    seq_info = "Removed " + seq_id + " " + records[seq_id].annotations["Species"]
                    seq_list.append(seq_info)
            print("")
            for seq in seq_list:
                print(seq)


def subset_records_with_regex(*header_terms, records, length=0, mode="exclude"):
    search_term = ""
    for term in header_terms:
        search_term += "\s" + term + "\d*-like|" + term + "\d*|"

    regex = re.compile(search_term[0:len(search_term) - 1])

    new_records = {}
    for record in records.values():
        if mode == "exclude":
            if len(str(record.seq)) > length and not regex.search(record.description):
                new_records[record.name] = record
        elif mode == "include":
            if len(str(record.seq)) > length and regex.search(record.description):
                new_records[record.name] = record

    return new_records


def exclude_character(records, character, mode="exclude"):
    new_records = {}
    for record in records.values():
        if mode == "exclude":
            if character not in record.seq:
                new_records[record.name] = record
        elif mode == "include":
            if character in record.seq:
                new_records[record.name] = record
    return new_records


def build_species_count(*header_terms, records, length=0, mode="exclude"):
    species_count = defaultdict(list)
    for record in records.values():
        if mode == "exclude":
            if len(str(record.seq)) > length and not any(term in record.description for term in header_terms):
                if "[" in record.description:
                    name = record.description.split("[")[1].split("]")[0]
            else:
                print(record.description)
            species_count[name].append(record.id)
    return species_count


def count_ids(records, min_length=0):
    count = 0
    for record in records.values():
        if len(record) > min_length:
            count += len(record)
    return count


def get_species_names(records, min_length=0, counts=False):
    species_names = []
    for k,v in records.items():
        if len(v) > min_length:
            species_names.append(k) if not counts else species_names.append(k + " " + str(len(v)))
    return species_names


def plot_record_number(records, plot_type, min_length=0):
    plot_records = {}
    # If we need to restrict the records to plot based on a minimum number
    if min_length > 0:

        for k, v in records.items():
            if len(v) >= min_length:
                plot_records[k] = v
    else:
        plot_records = records

    if plot_type == "Bar":

        data = [Bar(x=list(plot_records.keys()),
                    y=[len(x) for x in plot_records.values()])]

    elif plot_type == "Pie":
        trace = go.Pie(labels=plot_records.keys(), values=plot_records.values())

        py.iplot([trace], filename='basic_pie_chart')

    else:
        print("Invalid plot type selected. Choose 'Bar' or 'Pie'")

    return data


def map_dict_to_records(records, full_dict={}, unique=False):
    # If full_dict hasn't been provided, take records as the full dictionary
    if not full_dict:
        full_dict = records
    out_records = []
    for record in records.values():
        out_records.append(full_dict[record.id])
    return out_records


def map_species_dict_to_records(records, full_dict, unique=False):
    out_records = []
    for id_list in records.values():
        if unique:
            out_records.append(full_dict[id_list[0]])
        else:
            for seq_id in id_list:
                out_records.append(full_dict[seq_id])
    return out_records


def map_list_to_records(records, full_dict, unique=False):
    out_records = []
    for seq_id in records:
        if unique and seq_id in full_dict:
            out_records.append(full_dict[seq_id])
        elif seq_id in full_dict:
            out_records.append(full_dict[seq_id])
    return out_records


def write_fasta(records, filename):
    with open(filename, "w") as handle:
            SeqIO.write(records, handle, "fasta")


def compare_fasta(records1, records2):
    for x in records1.keys():
        if x not in records2.keys():
            print (x)


def replace_words(changes, file):
    for k, v in changes:
        file.replace(k, v)
    return file


def get_different_ids(records1, records2):
    ids1 = set(records1.keys())
    ids2 = set(records2.keys())
    return ids1 - ids2


def get_different_records(records1, records2):
    ids = get_different_ids(records1, records2)
    different_dict = {k: records1[k] for k in ids}
    return different_dict


def subset_on_motif(records, motif, retain_motif_seqs = True):
    """
    Subset a set of records based on the presence of a motif.

    :param motif: The motif to check for
    :param exclude: Whether we should retain sequences with the motif or exclude them
    :return: The subseted records based on the presence of the motif
    """
    subset_records = {}

    pattern=re.compile(motif)
    count = 1

    for record in records.values():
        if (len(pattern.findall(str(record.seq))) > 0 and retain_motif_seqs) or (len(pattern.findall(str(record.seq))) == 0 and not retain_motif_seqs):
            subset_records[record.id] = record
        else:
            print(count)
            print(record.id, pattern.findall(str(record.seq)))
            count +=1

    return subset_records

def correct_phylip_tree(records_path, phylip_file_path, outpath):
    """
    Using a dictionary mapping full id to PHYLIP id, change a PHYLIP annotated tree back to the full ID
    :param records: The alignment we will retrieve the original IDs from
    :param phylip_file: The tree with the truncated PHYLIP ids
    :return:
    """
    records = utilities.load_sequences(records_path)
    phylip_tree = utilities.load_tree(phylip_file_path)
    phylip_correction_dict = generate_phylip_correction_dictionary(records)
    print (phylip_correction_dict)
    # Temporary correction to non-unique names
    phylip_correction_dict["LIONS"] = "XP_021429049.1"
    phylip_correction_dict["CLAMS"] = "XP_018919738.1"

    for node in phylip_tree:
        print (node.name)
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
        if record.name[0:10] in phylip_correction_dict and not outpath:
            print("The PHYLIP names are non-unique")
            print (record.name[0:10])
            print (record.name)
        elif outpath:
            pass
        else:
            phylip_correction_dict[record.name[0:10]] = record.name

    return phylip_correction_dict



# correct_phylip_tree("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and "
#                     "bacteria/180312_fifty_percent_identity/Complete_sequences/PHYLIP_files/"
#                     "2U1_complete_probcons_removed.aln", "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding "
#                                                       "plants fungi nematodes insects and bacteria/"
#                                                       "180312_fifty_percent_identity/Complete_sequences/PHYLIP_files/"
#                                                       "PhyML_ProbCons.nwk", "/Users/gabefoley/Dropbox/PhD/Projects/2U1/"
#                                                                          "2U1_2018/Excluding plants fungi nematodes "
#                                                                          "insects and bacteria/180312_fifty_percent"
#                                                                          "_identity/Complete_sequences/"
#                                                                          "PHYLIP_files/PhyML_ProbCons_corrected.nwk")

# correct_phylip_tree("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and "
#                     "bacteria/180312_fifty_percent_identity/Complete_sequences/PHYLIP_files/"
#                     "2U1_complete_probcons_removed.aln", "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding "
#                                                       "plants fungi nematodes insects and bacteria/"
#                                                       "180312_fifty_percent_identity/Complete_sequences/PHYLIP_files/"
#                                                       "PhyML_ProbCons.nwk", "/Users/gabefoley/Dropbox/PhD/Projects/2U1/"
#                                                                          "2U1_2018/Excluding plants fungi nematodes "
#                                                                          "insects and bacteria/180312_fifty_percent"
#                                                                          "_identity/Complete_sequences/"
#                                                                          "PHYLIP_files/PhyML_ProbCons_corrected.nwk")


# full_record = SeqIO.to_dict(SeqIO.parse("files/candidates/regextest.fasta", "fasta"))
# full_record = SeqIO.to_dict(SeqIO.parse("/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes insects and bacteria/180312_fifty_percent_identity/Alignment_exclusion/2U1_50_percent_tagged.fasta",
# "fasta"))
# # motifSeqs = subset_on_motif(full_record, '[FW][SGNH].[GD][^F][RKHPT][^P]C[LIVMFAP][GAD]', retain_motif_seqs=True)
#
# original_motifSeqs = subset_on_motif(full_record, '[FW][SGNH].[GD][^F][RKHPT][^P]C[LIVMFAP][GAD]', retain_motif_seqs=False)
# motifSeqs = subset_on_motif(full_record, '[F]..[G]...C.[G]', retain_motif_seqs=False)
# # FXXGXbXXCXG
#
# print (len(original_motifSeqs))
# print (len(motifSeqs))
# for x in original_motifSeqs:
#     if x not in motifSeqs:
#         print (x)
# for record in motifSeqs.values():
# #     print(record.name)
#
# write_fasta(motifSeqs.values(), "/Users/gabefoley/Dropbox/PhD/Projects/2U1/2U1_2018/Excluding plants fungi nematodes
# insects and bacteria/180312_fifty_percent_identity/Alignment_exclusion/2U1_50_percent_motif.fasta" )

# seq1 = "NPS"
# seq2 = "NASNAS"
# seq3 = "NAT"
#
# pattern = re.compile('N[^P][ST]')
#
# if len(pattern.findall(seq2)) > 0:
#     print('found some')


# only_2U1_records = subset_records("2U1", "2U1-like", records=full_record, length=400, mode="include")
#
# print(only_2U1_records)

# test_records = subset_records_with_regex("2D", "2B", records=full_record, mode="include")
# for item in full_record:
#     print(full_record[item].description)
#
# print("---------------------------------")
#
# for item in test_records:
#     print(test_records[item].description)
