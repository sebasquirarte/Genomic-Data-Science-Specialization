# Comparison betwreen Naive and Boyer–Moore string-search algorithms
# Part of Genomic Data Science Specialization - Algorithms for DNA Sequencing by Johns Hopkins University through Coursera
# Sebastian Quirarte | sebastianquirajus@gmail.com | 29 Oct 22

### FUNCTIONS ###
# Parses a DNA r eference genome from a file in the FASTA format
def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

# Naive exact matching algorithm (without reverse complement)
def naive(p, t):
    occurrences = []
    alignments = 0
    comparisons = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        alignments += 1
        match = True
        for j in range(len(p)):  # loop over characters
            comparisons += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, alignments, comparisons

# Naive matching algorithm (without reverse complement, allows for n mismatches)
def naive_mm(p, t, n):  # pattern, text, allowed mismatches
    occurrences = []
    alignments_count = 0
    comparisons_count = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mismatches = 0
        match = True
        alignments_count += 1
        for j in range(len(p)):  # loop over characters
            comparisons_count += 1
            if t[i+j] != p[j]:  # compare characters
                mismatches += 1
                if mismatches > n:
                    match = False
                continue
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, alignments_count, comparisons_count

# Boyer-Moore preprocessing
__author__ = "Ben Langmead"
import unittest

def z_array(s):
    """ Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)

    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break
    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1
    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]| """
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """ Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. """
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """ Given a full match of P to T, return amount to shift as
        determined by good suffix rule. """
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab

class BoyerMoore(object):
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, p, alphabet='ACGT'):
        # Create map from alphabet characters to integers
        self.amap = {alphabet[i]: i for i in range(len(alphabet))}
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule at offset i """
        assert c in self.amap
        assert i < len(self.bad_char)
        ci = self.amap[c]
        return i - (self.bad_char[i][ci]-1)

    def good_suffix_rule(self, i):
        """ Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule. """
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        """ Return amount to shift in case where P matches T """
        return len(self.small_l_prime) - self.small_l_prime[1]

# Boyer Moore matching algorithm
def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    alignments = 0
    comparisons = 0
    while i < len(t) - len(p) + 1:
        alignments += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            comparisons += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, alignments, comparisons

# Implementating a K-mer Index

__author__ = "Ben Langmead"

import bisect

class Index(object):
    """ Holds a substring index for a text T """
    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k # k-mer length (k)
        self.index = []
        for i in range(len(t)-k+1): # for each k-mer
            self.index.append((t[i:i+k], i)) # add (k-mer, offset) pair
        self.index.sort()

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k] # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1)) # binary search
        hits = []
        while i < len(self.index): # # collect matching index entries
            if self.index[i][0] != kmer:
                break # end of multimap equal range
            hits.append(self.index[i][1])
            i += 1
        return hits

# Partial matching algroithm implementation using an exact matching algorithm (Indexing)
# Combined with the pigeon hole principle to allow up to k mismatches of pattern in text
def approximate_match_index(p, t, k): #pattern, text, mismatches
    segment_length = round(len(p) // (k+1))
    hits = 0
    all_matches = set()
    index = Index(t, 8) # built on 8-mers
    for i in range(k+1):
        start = i * segment_length
        end = min((i+1) * segment_length, len(p))
        matches = index.query(p[start:end])
        hits += len(matches)
        for m in matches:
            text_offset = m - start
            if text_offset < 0 or (text_offset + len(p)) > len(t):
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[text_offset + j]:
                    mismatches += 1
                    if mismatches > k:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[text_offset + j]:
                    mismatches += 1
                    if mismatches > k:
                        break
            if mismatches <= k:
                all_matches.add(text_offset)
    return list(all_matches), hits

### CODE ###
# Sequence from fasta file, pattern, and allowed mismatches
seq = readGenome('chr1.GRCh38.excerpt.fasta')
p = 'GGCGCGGTGGCTCACGCCTGTAAT'
n = 2

# Variables
occurrences_naiveE, alignments_naiveE, comparisons_naiveE = naive(p, seq)
occurrences_naivemm, alignments_naivemm, comparisons_naivemm = naive_mm(p, seq, n)
occurrences_BM, alignments_BM, comparisons_BM = boyer_moore(p, BoyerMoore(p), seq)
occurrences_AM, hits_AM = approximate_match_index(p, seq, n)

# Naive exact search
print("\nNAIVE EXACT SEARCH:")
print("match count:",str(len(occurrences_naiveE)),"\nalignments:",alignments_naiveE,"\ncomparisons:",comparisons_naiveE)
while True:
    print_indexes = input("See all indexes? (y/n): ")
    if print_indexes == "y" or print_indexes == "Y":
        print(occurrences_naiveE)
        break
    elif print_indexes == "n" or print_indexes == "N":
        break
    else:
        print("Invalid input.")

# Naive w/ allowed mismatches
print("\nNAIVE SEARCH:")
print("match count:",str(len(occurrences_naivemm)),"\nalignments:",alignments_naivemm,"\ncomparisons:",comparisons_naivemm)
while True:
    print_indexes = input("See all indexes? (y/n): ")
    if print_indexes == "y" or print_indexes == "Y":
        print(occurrences_naivemm)
        break
    elif print_indexes == "n" or print_indexes == "N":
        break
    else:
        print("Invalid input.")

# Bayer-Moore search
print("\nBOYER-MOORE SEARCH:")
print("match count:",str(len(occurrences_BM)),"\nalignments:",alignments_BM,"\ncomparisons:",comparisons_BM)
while True:
    print_indexes = input("See all indexes? (y/n): ")
    if print_indexes == "y" or print_indexes == "Y":
        print(list(occurrences_BM),"\n")
        break
    elif print_indexes == "n" or print_indexes == "N":
        break
    else:
        print("Invalid input.")

# Match Index Search
print("\nMATCH INDEX SEARCH:")
approximate_match_index(p, seq, n)
print("match count:",str(len(occurrences_AM)),"\nindexes:",hits_AM)
while True:
    print_indexes = input("See all indexes? (y/n): ")
    if print_indexes == "y" or print_indexes == "Y":
        print(list(occurrences_AM),"\n")
        break
    elif print_indexes == "n" or print_indexes == "N":
        break
    else:
        print("Invalid input.")



##### TODO 1. 2. Check what a class is, 3. Check what __init__ is and what main() is
