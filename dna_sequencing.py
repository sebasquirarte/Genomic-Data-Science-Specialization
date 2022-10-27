# Multiple functions related to DNA sequencing
# Part of the 'Genomic Data Science' Specialization | Course 3: Algorithms for DNA Sequencing
# Sebastian Quirarte | sebastianquirajus@gmail.com | 25 Oct 22

# Takes a DNA string and returns its reverse complement
def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

# Parses a DNA reference genome from a file in the FASTA format
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

# Parses the read and quality strings from a FASTQ file containing sequencing reads
def readFastq(filename):
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

# Turns	Q into Phred+33	ASCII-encoded quality
def QtoPhred33(Q):
		return chr(Q	+ 33)

# Turns Phred+33 ASCII-encoded quality into Q
def phred33ToQ(qual):
		return ord(qual)-33

# Naive exact matching algorithm (without reverse complement)
def naive(p, t):  # pattern, text
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    print(len(occurrences),"matches found")
    print("Leftmost match index:",min(occurrences))
    print("Rightmost match index:",max(occurrences))
    while True:
        print_indexes = input("\nSee all indexes? (y/n): ")
        if print_indexes == "y" or print_indexes == "Y":
            print(occurrences)
            break
        elif print_indexes == "n" or print_indexes == "N":
            break
        else:
            print("Invalid input.")

# Naive exact matching algorithm (without reverse complement)
def naive_with_rc(p, t):  # pattern, text
    occurrences = []
    rcomplement_p = reverseComplement(p)
    print("\nReverse Complement: ",rcomplement_p)
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    if p != rcomplement_p:
        print("Reverse complement IS NOT the same")
        for i in range(len(t) - len(rcomplement_p) + 1):  # loop over alignments
            match = True
            for j in range(len(rcomplement_p)):  # loop over characters
                if t[i+j] != rcomplement_p[j]:  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)  # all chars matched; record
    else:
        print("Reverse complement IS the same")
    print(len(occurrences),"matches found")
    print("Leftmost match index:",min(occurrences))
    print("Rightmost match index:",max(occurrences))
    while True:
        print_indexes = input("\nSee all indexes? (y/n): ")
        if print_indexes == "y" or print_indexes == "Y":
            print(occurrences)
            break
        elif print_indexes == "n" or print_indexes == "N":
            break
        else:
            print("Invalid input.")

# Naive exact matching algorithm (without reverse complement, allows up to 2 mismatches)
def naive_2mm(p, t):  # pattern, text
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatches = 0
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatches += 1
                if mismatches > 2:
                    match = False
                continue
        if match:
            occurrences.append(i)  # all chars matched; record
    print("\n" + str(len(occurrences)),"matches found")
    print("Leftmost match index:",min(occurrences))
    print("Rightmost match index:",max(occurrences))
    while True:
        print_indexes = input("\nSee all indexes? (y/n): ")
        if print_indexes == "y" or print_indexes == "Y":
            print(occurrences)
            break
        elif print_indexes == "n" or print_indexes == "N":
            break
        else:
            print("Invalid input.")

# Creates histogram of qualities and their frequency
def createHist(qualities):
    hist = [0] * 50
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] += 1
    return hist

# Opens file
while True:
    file = input("Enter file name: " )
    if file == "":
        break
    try:
        fileh = open(file)
        break
    except:
        print("No such file in directory.")

# For .fa files
pattern = input("Enter DNA pattern (ALL CAPS): ")
genome = readGenome(file)
#naive(pattern, genome)
naive_with_rc(pattern, genome)
#naive_2mm(pattern, genome)

# For .fastaq files
# sequences, qualities = readFastq(file)

# Plots histogram
#import matplotlib.pyplot as plt
#h = createHist(qualities)
#plt.bar(range(len(h)),h)
#plt.show()
