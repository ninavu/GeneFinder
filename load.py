import readline

# a function that takes the name of the FASTA file and returns all the DNA 
def loadSeq(fileName):
    fin = open(fileName, "r")
    first_line = fin.readline() # skip the first line
    lines = fin.readlines()
    
    str = ""
    for line in range(len(lines)):
        str += lines[line].strip()
    fin.close()
    return str

# def main():
#     # print(loadSeq("U81861.fa"))
#     print(loadSeq("X73525.fa"))

# main()