from Bio import SeqIO, Entrez


# handle = SeqIO.parse("/Users/gabefoley/Dropbox/Code/Python Workspace/digital_notebook/files/full.fasta", "fasta")
handle = SeqIO.parse("/Users/gabefoley/Dropbox/YenA2_psiblast_top100_seqdump.fasta", "fasta")

species_set = set()
species_list = []
unique_strain_list = []
unique_strain_set = set()
for seq in handle:
    name = seq.description.split("[")[1].split("]")[0]
    species_set.add(name)
    if name not in species_list:
        species_list.append(name)
    if name not in unique_strain_list and "sp." not in name:
        unique_strain_list.append(name)

        # if len(name.split()) > 1:
        #     if name.split()[0] not in unique_strain_set:
        #         unique_strain_list.append(name)
        #         unique_strain_set.add(name.split()[0])
idx = 1

print ("These are the hits for unique strains")
for name in unique_strain_list:
    if idx <= 102:
        print (name)
        idx += 1
        if len(name.split(" ")) == 1:
            print ("***")
print ("\n**************************************")
print ("These are the hits for unique species")
idx = 1
for name in species_list:
    if idx <= 102:
        print(name)
        idx += 1
        if len(name.split(" ")) == 1:
            print("***")
# for name in species_set:
#     if idx <= 102:
#         print (name)
#         idx += 1
#         if len(name.split(" ")) == 1:
#             print ("***")

# while species_set:
#     with open("/Users/gabefoley/Dropbox/Code/Python Workspace/digital_notebook/files/species_" + str(idx) + ".txt", "a") as file:
#         for count in range (0,50):
#             if (species_set):
#                 file.write(species_set.pop() + "\n")
#         idx +=1
