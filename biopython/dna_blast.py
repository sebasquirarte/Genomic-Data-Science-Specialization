# Runs blast (NCBI) on dna sequence
# Sebastian Quirarte | sebastianquirajus@gmail.com
# 24 Oct 22

from Bio.Blast import NCBIWWW, NCBIXML

# Insert sequence here
result_handle = NCBIWWW.qblast("blastn", "nt", "tgggcctcatatttatcctatataccatgttcgtatggtggcgcgatgttctacgtgaatccacgttcgaaggacatcataccaaagtcgtacaattaggacctcgatatggttttattctgtttatcgtatcggaggttatgttcttttttgctctttttcgggcttcttctcattcttctttggcacctacggtagag")
blast_record = NCBIXML.read(result_handle)

# Modify e value if needed
e_value_thresh = 0.01
count = 0
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < e_value_thresh:
            count += 1
            print("\n---------Alignment---------")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)

print("\nThere are",count,"similar sequences in Blast output")
