from Bio import SeqIO
from Bio.Restriction import NotI


record = SeqIO.read("sequence.fasta", "fasta")
print(record)
sites = NotI.search(record.seq)
print(sites, len(sites))
