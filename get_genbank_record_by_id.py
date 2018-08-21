import os, sys
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "carolina.grandellis@tsl.ac.uk"

idlistfile=open(sys.argv[1])
ids=[]
for line in idlistfile:
	line=line.rstrip()
	if line == "": continue
	ids.append(line)

idlistfile.close()
a=0
while a <=len(ids) + 100:
	list_of_ids=",".join(ids[a:a+100])
	#print(list_of_ids)
	with Entrez.efetch(db="protein", rettype="gb", retmode="text", id=list_of_ids) as handle:
		for record in SeqIO.parse(handle, "gb"):
			for (index, feature) in enumerate(record.features):
				if feature.type == 'CDS':
					geneid=feature.qualifiers['db_xref'][0].replace('GeneID:', '')

					print(record.id + "\t" + feature.qualifiers['db_xref'][0] + "\thttps://www.ncbi.nlm.nih.gov/gene/" + geneid)
	a=a+100

exit(0)
