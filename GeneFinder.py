from tracemalloc import start
from aminoAcids import *
from load import *
from dna import *
import random

startCodon = codons[3]
stopCodons = codons[10]

# a function that collects all possible open reading frames in a DNA sequence
def oneFrame(DNA):
    frameList = []
    for i in range(0, len(DNA), 3):
        if DNA[i:i+3] in startCodon: # only start adding when start codon is found
            frame = ""
            j = i
            while (DNA[j:j+3] not in stopCodons and j < len(DNA)):
                frame += DNA[j:j+3]
                j += 3
            frameList.append(frame)
    return frameList

# a function that chooses the longest ORF among the a set of nested ORFs
def oneFrameV2(DNA):
    frameList = []
    i = 0
    while (i < len(DNA)):
        if (DNA[i:i+3] not in startCodon):
            i += 3
        else:
            j = i
            frame = ""
            while (DNA[j:j+3] not in stopCodons and j < len(DNA)):
                frame += DNA[j:j+3]
                j += 3
            frameList.append(frame)
            i = j   
    return frameList 

# a function that finds the longest ORF of all 3 reading frames
def longestORF(DNA):
    longORF = ""
    for i in range(3):
        ORF = oneFrameV2(DNA[i:])
        for j in range (len(ORF)):
            if len(ORF[j]) > len(longORF):
                longORF = ORF[j]
    return longORF

# a function that finds the longest ORF between the DNA stand and its reverse complement
def longestORFBothStrands(DNA):
    ORF_DNA = longestORF(DNA) 
    ORF_rcDNA = longestORF(reverseComplement(DNA))
    if len(ORF_DNA) >= len(ORF_rcDNA):
        return ORF_DNA
    else:
        return ORF_rcDNA    

# a function that takes a list as an input and returns back the string
def collapse(L):
    str = ""
    for i in range(len(L)):
        str += L[i]
    return str

# a function that makes garbage sequences, finds the longest ORF, and returns its length
def longestORFNoncoding(DNA, numReps):
    listDNA = list(DNA)
    maxORF = ""
    for i in range(numReps):
        random.shuffle(listDNA)
        newDNA = collapse(listDNA)
        ORF = longestORFBothStrands(newDNA)
        if len(ORF) >= len(maxORF):
            maxORF = ORF
    return (len(maxORF))

# a function that identifies all the ORFs in the unshuffled DNA and returns them as a list
def findORFs(DNA):
    listORF = []
    for i in range(3):
        ORF = oneFrameV2(DNA[i:])
        for j in range (len(ORF)):
            listORF.append(ORF[j])
    return listORF

# a function that searches both the forward and reverse complement strands for ORFs
# and returns a list with all the ORFs found 
def findORFsBothStrands(DNA):
    allORFs = []
    ORF_DNA = findORFs(DNA)
    ORF_rcDNA = findORFs(reverseComplement(DNA))

    for i in range(len(ORF_DNA)):
        allORFs.append(ORF_DNA[i])

    for i in range(len(ORF_rcDNA)):
        allORFs.append(ORF_rcDNA[i])
    return allORFs

# a function that returns the beginning and end coordinates of an ORF in DNA
def getCoordinates(orf, DNA):
    listCoordinates = []
    rcORF = reverseComplement(orf)

    # ORF falls on the forward strand
    startC = DNA.find(orf)
    if startC == -1:
        startC = DNA.find(rcORF)
    listCoordinates.append(startC)

    endC = startC + len(orf)
    listCoordinates.append(endC)
    return listCoordinates

# a function that identifies ORFs longer than minLen and returns a list with info of each one
def geneFinder(DNA, minLen):
    listORF = findORFsBothStrands(DNA)
    finalOutputList = []

    for i in range(len(listORF)):
        curORF = listORF[i]
        if len(curORF) > minLen:
            infoORF = [getCoordinates(curORF, DNA)[0],
                       getCoordinates(curORF, DNA)[1],
                       codingStrandToAA(curORF)]
            finalOutputList.append(infoORF)
    finalOutputList.sort()    
    return finalOutputList

def printGenes(geneList):
    print("*** Output of geneFinder ***")
    for i in range(len(geneList)):
        print("---Gene number: " + str(i+1))
        print("Start Coordinate: " + str(geneList[i][0]))
        print("End Coordinate: " + str(geneList[i][1]))
        print("Corresponding amino acids: " + geneList[i][2])
    return 

def main():
    # ---Test functions
    # print(oneFrame("CCCATGTTTTGAAAAATGCCCGGGTAAA"))
    # print(oneFrame("CCATGTTAGAAATGCCC"))
    # print(oneFrame("ATGCCCATGGGGAAATTTTGACCC"))

    # print(oneFrameV2("ATGCCCATGGGGAAATTTTGACCC"))
    # print(oneFrameV2("TGCCCTAACATGAAAATGACTTAGG"))

    # print(longestORF("ATGAAATAG"))
    # print(longestORF("CATGAATAGGCCCA"))
    # print(longestORF("CTGTAA"))
    # print(longestORF("ATGCCCTAACATGAAAATGACTTAGG"))

    # print(longestORFBothStrands("CTATTTCATG"))

    # X73525 = loadSeq("X73525.fa")
    # print(longestORFNoncoding(X73525, 50))
    # print(longestORFNoncoding(X73525, 50))

    # print(findORFs("ATGGGATGAATTAACCATGCCCTAA"))
    # print(findORFs("GGAGTAAGGGGG"))

    # print(findORFsBothStrands("ATGAAACAT"))

    # print(getCoordinates("GTT", "ACGTTTCGA"))
    # print(getCoordinates("CGAA", "ACGTTCGA"))

    # print(geneFinder("ATGAAATAG", 2))

    # Gene Finder
    X73525 = loadSeq("X73525.fa")
    minLen = longestORFNoncoding(X73525, 1500)
    print("Threshold value / minLen: " + str(minLen))
    geneList = geneFinder(X73525, minLen)
    printGenes(geneList)

main()