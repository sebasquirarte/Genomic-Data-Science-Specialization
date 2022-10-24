# Write a Python program that takes as input a file containing DNA sequences in multi-FASTA format,
# and computes the answers to the questions in 'instructions.txt'.
# Sebastian Quirarte | sebastianquirajus@gmail.com | 22 Oct 22 

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

# Creates dictionay and stores values
sequences = {}
for line in fileh:
    line = line.rstrip()
    # Checks if line is header
    if line.startswith(">"):
        words = line.split()
        name = words[0][1:]
        sequences[name] = ""
    # Adds sequence if line is not a header
    else:
        sequences[name] = sequences[name] + line

sequences_length = {}

# Prints names, secuences, and/or length
for name, seq in sequences.items():
    # Print name, full sequences and length
    #print("\n" + name + "\n" + seq + "\nLength: " + str(len(seq)))
    # Print name and length
    #print("\n" + name + "\nLength: " + str(len(seq)))
    # Adds values to list of lengths
    sequences_length[name] = len(seq)

# Determines largest and smallest sequence
largest = max(sequences_length.values())
smallest = min(sequences_length.values())
largest_name = list(sequences_length.keys())[list(sequences_length.values()).index(largest)]
smallest_name = list(sequences_length.keys())[list(sequences_length.values()).index(smallest)]

# Prints values and divison line
print("\nSecuences: " + str(len(sequences)))
print("Largest:",largest_name,"Count:",largest)
print("Smallest:",smallest_name,"Count:",smallest,"\n\n______________________________________ORFs______________________________________\n")

#ORFs functions
orfs1 = []
orfs2 = []
orfs3 = []
all_orfs = []

def find_orf(sequence, n):
    print("--READING FRAME " + str(n) + "--")
    # Finds start and stop codon indexes
    start_index = []
    stop_index = []
    for i in range(n-1, len(sequence), 3):
        if sequence[i:i+3] == "ATG":
            start_index.append(i)
        if sequence[i:i+3] in ["TAA", "TAG", "TGA"]:
            stop_index.append(i+3)
    # Finds ORFs
    orfs = []
    position = []
    marker = 0
    for i in range(0, len(start_index)):
        for j in range (0, len(stop_index)):
            if start_index[i] < stop_index[j] and start_index[i] > marker:
                orfs.append(sequence[start_index[i]:stop_index[j]])
                all_orfs.append(sequence[start_index[i]:stop_index[j]])
                position.append(start_index[i])
                mark = stop_index[j] + 3
                if n == 1:
                    orfs1.append(sequence[start_index[i]:stop_index[j]])
                elif n == 2:
                    orfs2.append(sequence[start_index[i]:stop_index[j]])
                elif n == 3:
                    orfs3.append(sequence[start_index[i]:stop_index[j]])
                break
    try:
        longest_index = (orfs.index(max(orfs, key=len)))
        shortest_index = (orfs.index(min(orfs, key=len)))
    except:
        pass
    orfs.sort(key=len)
    try:
        print("Longest ORF: " + str(len(orfs[-1])),"Index: " + str(position[longest_index] + 1))
    except:
        print("None")
    try:
        print("Shortest ORF: " + str(len(orfs[0])),"Index: " + str(position[smallest_index] + 1))
    except:
        print("None")

# Finds ORFs
for key, sequence in sequences.items():
    print(key)
    find_orf(sequence,1)
    find_orf(sequence,2)
    find_orf(sequence,3)

# Sorts ORFs
all_orfs.sort(key=len)
orfs1.sort(key=len)
orfs2.sort(key=len)
orfs3.sort(key=len)

# Prints longest and shortest from all ORFs and reading frames
print("-----ALL ORFS-----\n" + "Longest:",str(len(all_orfs[-1])),"\nShortest:",str(len(all_orfs[0])))
print("-----READING FRAME 1-----\n" + "Longest:",str(len(orfs1[-1])),"\nShortest:",str(len(orfs1[0])))
print("-----READING FRAME 2-----\n" + "Longest:",str(len(orfs2[-1])),"\nShortest:",str(len(orfs2[0])))
print("-----READING FRAME 3-----\n" + "Longest:",str(len(orfs3[-1])),"\nShortest:",str(len(orfs3[0])))

fileh.close()
