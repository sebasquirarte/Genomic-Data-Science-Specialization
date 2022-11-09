# Function that determines the edit distance between a sequence and a genome using a matrix
# Part of Genomic Data Science Specialization - Algorithms for DNA Sequencing by Johns Hopkins University through Coursera
# Sebastian Quirarte | sebastianquirajus@gmail.com | 9 Nov 22

# Opens and reads genome from 'fasta' file
def readFASTA(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

# Opens and reads genome from 'fastq' file
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

def editDistance(p, t): # pattern, text
    m, n = len(p), len(t)
    dp = [[0 for j in range(n+1)] for i in range(m+1)]
    for i in range(m+1): dp[i][0] = i # init first column by distance from empty string
    for i in range(1, m+1):
        for j in range(1, n+1):
            dp[i][j] = min(
                dp[i-1][j-1] + int(p[i-1] != t[j-1]),
                dp[i-1][j] + 1,
                dp[i][j-1] + 1,)
    return min(dp[m])

def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def overlap_all_pairs(reads, k, map={}): # k is the minimum overlap length
    def get_kmers(read, k):
        res = set()
        for i in range(0, len(read)-k+1):
            res.add(read[i:i+k])
        return res
    for read in reads:
        kmers = get_kmers(read, k)
        for kmer in kmers:
            if not kmer in map.keys():
                map[kmer] = set()
            map[kmer].add(read)
    pairs = []
    for head in reads:
        kmer = head[-k:]
        candidates = map[kmer]
        for tail in candidates:
            if (not head == tail and overlap(head, tail, k)):
                pairs.append((head, tail))
    return pairs

genome = readFASTA("chr1.GRCh38.excerpt.fasta")
p1 = "GCTGATCGATCGTACG"
p2 = "GATTTACCAGATTGAG"

print("Reading chr1.GRCh38.excerpt.fasta...")
print("Pattern:",p1,"\nEdit Distance:",(editDistance(p1, genome)))
print("\nPattern:",p2,"\nEdit Distance:",(editDistance(p2, genome)))

reads, _ = readFASTQ("ERR266411_1.for_asm.fastq")

print("\nReading ERR266411_1.for_asm.fastq...")
pairs = overlap_all_pairs(reads, 30)
print("Pairs (Edges in Graph):",str(len(pairs)),"\nNodes w outgoing edges:",str((len(set(pair[0] for pair in pairs)))))
