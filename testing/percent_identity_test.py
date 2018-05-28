import alignment
import utilities
import glob

working_dir = "/Users/gabefoley/Dropbox/Code/Python Workspace/digital_notebook/notebooks/Latest exon testing/All_Full_aligned/"
# working_dir = "/Users/gabefoley/Dropbox/Code/Python Workspace/digital_notebook/testing/example"
files = glob.glob(working_dir + "/*.fasta")

with open("./output.csv", "w+") as out_file:
    out_file.write("Family, Mean Percent Identity, Standard Deviation of Percent Identity\n")
    for file in files:
        name = file.split("/")[-1].split(".")[0]
        print ("Working on %s" % (name))

        alignment_file = utilities.load_alignment(file)

        percent_identities_array = alignment.get_percent_identity_of_alignment(alignment_file, count_gaps=False)

        print("Percent identity is ", percent_identities_array.mean())
        print("Standard deviation is ", percent_identities_array.std())
        out_file.write("%s,%s,%s\n" % (name, percent_identities_array.mean(), percent_identities_array.std()))


