from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt


seqs = SeqIO.parse('16S-UDb_sequences.fasta','fasta')
taxs = open('16S-UDb_MOTHUR_taxonomy.txt','r')
lines = taxs.readlines()

# make a matrix containing several list that contains each id and taxanomy
# ex) [[id1,kingdom,phylum,...,species], [id2,kingdom,phylum,...,species],...]
taxList = []
for line in lines:
    line = line.strip()
    id,semi_tax = line.split('\t')[0], line.split('\t')[1]
    taxs = semi_tax.split(';')
    id_tax_conc = []
    id_tax_conc.append(id)
    for tax in taxs:
        id_tax_conc.append(tax)
    taxList.append(id_tax_conc)

# dictionary list containing id and its sequence
seqDict = {}
for seq_record in seqs:
    seqDict[seq_record.id] = str(seq_record.seq[:5])

# add sequence into the matrix
for tax in taxList:
    tax.append(seqDict[tax[0]])

# delete errors in species' name
result = []
for tax in taxList:
    if len(tax) == 9:
        result.append(tax)

# dataframing
info_column = ['ID', 'Kingdom','Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'seq']
info_df = pd.DataFrame(data=result, columns=info_column)
print(info_df.head(100))


# visualization of species distribution
fam_amount = dict()
for fam in info_df["Species"]:
    if fam in fam_amount:
        fam_amount[fam] += 1
    else:
        fam_amount[fam] = 1

myList = fam_amount.items()
myList = sorted(myList)
x, y = zip(*myList)

plt.plot(x, y)
# plt.show()
species = pd.DataFrame(myList, columns=["Species", "Amount"])



