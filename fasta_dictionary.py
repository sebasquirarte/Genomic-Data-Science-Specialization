# fasta_dictionary - Reads a 'fasta' file and creates a dictionary containing all the DNA or AA's sequences in the file
# 22 Oct 22 | Sebastian Quirarte | sebastianquirajus@gmail.com

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

# Prints stored values
for name, seq in sequences.items():
    print("\n" + name + "\n" + seq)

print("")
fileh.close()
