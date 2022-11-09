# Functions that determine the overlaps and Shortest Common String (SCS) to assemble a genome from a set of sequences
# Part of Genomic Data Science Specialization - Algorithms for DNA Sequencing by Johns Hopkins University through Coursera
# Sebastian Quirarte | sebastianquirajus@gmail.com | 9 Nov 22

import itertools

# Opens and reads 'fastq' file
def readFASTQ(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def scss(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    shortest_cnt = 1
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
            shortest_cnt = 1
        elif len(sup) == len(shortest_sup):
            shortest_cnt += 1
    return shortest_sup, shortest_cnt  # return shortest

def pick_max_overlap(reads, k):
    best_a, best_b, best_len = None, None, 0
    for a, b in itertools.permutations(reads, 2):
        length = overlap(a, b, min_length=k)
        if length > best_len:
            best_a, best_b, best_len = a, b, length
    return best_a, best_b, best_len

def greedy_scss(reads, k):
    while True:
        a, b, olen = pick_max_overlap(reads, k)
        if olen == 0:
            break
        reads.remove(a)
        reads.remove(b)
        reads.append(a + b[olen:])
    return ''.join(reads) # append all non-overlaps onto eachother and return the concatenated string

# shortest_common, frequency = scss(['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'])
# print(len(shortest_common))
# print(frequency)

reads, _ = readFASTQ('mysteryvirus.fastq')
print("Reading file...")
genome = greedy_scss(reads, 30)
print("Genome:",genome)
print("Length of genome:",str(len(genome)))
print("A count:",str(mystery_genome.count('A')))
print("B count:",str(mystery_genome.count('T')))
