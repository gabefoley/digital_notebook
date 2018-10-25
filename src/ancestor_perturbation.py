import sys
import fasta
import os
import utilities
import subprocess
import alignment
import time

devnull = open(os.devnull, 'w')

# CYP2U1_CYP2R1_Realigned_PPG.aln CYP2U1_CYP2R1_Realigned_Met_Rerooted.nwk


# Grab all the keyword arguments

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: ancestor_perturbation <alignment> <tree> where', \
              "\n\t<alignment> is the multiple sequence alignment", \
              "\n\t<tree> is the phylogentic tree")
        sys.exit(1)
    alignment_file = utilities.load_sequences( os.getcwd() + "/" + sys.argv[1])
    tree = utilities.load_tree(os.getcwd() + "/" + sys.argv[2])



    # Create variables with just the names, not the extensions of the inputs or the input objects
    alignment_name = sys.argv[1].split(".")[0]
    tree_name = sys.argv[2].split(".")[0]

    align = '-align' in sys.argv[3:]
    reinfer = '-reinfer' in sys.argv[3:]




# Remove each sequence in turn
count = 1
start_time = time.time()

for seq in alignment_file.values():
    if seq.name == "XP_008105426.1":
        segment_start = time.time()
        print ("Working on sequence %d of %d" % (count, len(alignment_file.values())))

        # Remove sequence from FASTA
        filtered_aln = fasta.subset_records(seq.name, records=alignment_file)
        new_sequences_name = os.getcwd() + "/" + alignment_name + "_without_" + seq.name + ".fasta"
        new_alignment_name = os.getcwd() + "/" + alignment_name + "_without_" + seq.name + ".aln"
        new_tree_name = tree_name + "_without_" + seq.name + ".nwk"



        # If realign, realign the sequences
        if align:
            print ("Realigning")
            for seq_name in filtered_aln:
                filtered_aln[seq_name].seq = filtered_aln[seq_name].seq.ungap("-")


            fasta.write_fasta(filtered_aln.values(), new_sequences_name, True)

            # Create a new alignment
            new_alignment = alignment.align_with_mafft(new_sequences_name, localpair=True)
            # Save the alignment file
            fasta.write_fasta(new_alignment, new_alignment_name, True)
        # Else remove any gap only columns

        else:
            fasta.write_fasta(filtered_aln.values(), new_alignment_name ,True)
            subprocess.call(["trimal", "-in", new_alignment_name, "-out", new_alignment_name, "-noallgaps"])




        # If infer, infer a new tree

        if reinfer:
            print ("Inferring a new tree")
            tree_path = os.getcwd() + "/" + "RAxML_bestTree." + tree_name + "_without_" + seq.name + ".nwk"
            subprocess.call(["raxml", "-m", "PROTGAMMAIJTT", "-p", "23456", "-n", new_tree_name, "-s", new_alignment_name], stdout=devnull)

            while not os.path.exists(tree_path):
                time.sleep(1)

            if os.path.isfile(tree_path):


                pared_tree = utilities.load_tree(tree_path)
                pared_tree.set_outgroup("XP_020383892.1")
                pared_tree.write(outfile=new_tree_name)

                # print ("RAxML_bestTree." + tree_name + "_without_" + seq.name + ".nwk")
                # subprocess.call(["mv", tree_path, new_tree_name])




        # # Else remove sequence from tree

        else:
            subprocess.call(["java", "-jar", "Paretree.jar", "-del", seq.name, "-f", tree_name + ".nwk"], stdout=devnull)

            pared_tree = utilities.load_tree(os.getcwd() + "/" + tree_name + "_pared.nwk")
            pared_tree.write(outfile=new_tree_name)

            # subprocess.call(["mv",tree_name + "_pared.nwk", new_tree_name])

        output_folder = "output/" + "output_without_" + seq.name

        subprocess.call(["mkdir", output_folder], stdout=devnull, stderr=devnull)




        # # Input new alignment and new tree into bnkit

        print ("Inferring new ancestor")
        subprocess.call(["java", "-jar", "bnkit.jar", new_tree_name, "-s", new_alignment_name, "-mp", "-model", "JTT", "-inf", "Joint", "-p", "3", "-o", output_folder])

        while not os.path.exists(output_folder + "/_aln_full.fa"):
            time.sleep(1)

        if os.path.isfile(output_folder + "/_aln_full.fa"):

            ancestor_aln = utilities.load_sequences(output_folder + "/_aln_full.fa")

            filtered_ancestor = fasta.subset_records("N1_1", records=ancestor_aln, mode="include")

            with open(os.getcwd() + "/output/final_output.fasta", "a") as output_file:
                output_file.write(">" + new_alignment_name +"\n")
                output_file.write(str(filtered_ancestor["N1_1"].seq) + "\n")

            count +=1

            print (new_alignment_name)
            print (filtered_ancestor["N1_1"].seq)
            segment_end = time.time()

        print ("Finished sequence %s in %s" % (count, segment_end - segment_start))


print ("Finished. Total time was %s " % (time.time() - start_time))




