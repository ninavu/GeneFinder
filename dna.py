from aminoAcids import *

# a function that takes a DNA string as input
# and returns the sequence of its complementary DNA strand
def reverseComplement(DNA):
    # reverse the string
    DNA = DNA[::-1]

    # match with the the complementary strand by modifying a list
    l = list(DNA)
    for i in range(len(DNA)):
        if l[i] == 'A':
            l[i] = 'T'
        elif l[i] == 'T':
            l[i] = 'A'
        elif l[i] == 'G':
            l[i] = 'C'
        elif l[i] == 'C':
           l[i] = 'G'
    return(''.join(l))

# a function that takes a DNA sequence from the coding strand 
# and returns the corresponding amino acids as a string
def codingStrandToAA(DNA):
    if len(DNA) % 3 != 0:
        print("DNA sequence not divisible by 3!")
        return None
    else:
        # translate RNA to protein sequence by checking first - second - third base
        # using U in RNA as T in DNA
        AA = ""
        for i in range(0,len(DNA),3):
            bases = DNA[i:i+3]
            for j in range(len(codons)):
                for k in range(len(codons[j])):
                    if bases == codons[j][k]:
                        AA += aa[j]
        return AA

# def main():
#     # Test 1st function
#     print((reverseComplement("TTGAC"))) # GTCAA

#     # Test 2nd function
#     print((codingStrandToAA("AGTCCCGGGTTT")))  # SPGF
#     print((codingStrandToAA("ATGCAACAGCTC")))  # MQQL
#     print((codingStrandToAA("ATGCAACAGCTCT"))) # error

# main()