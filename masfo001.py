#-------------------------------------------------------------------------------
# Author:      amir
# Created:     10/25/2015
#
# Instructions:
#
# 1) Make sure to rename the file (studentNetId.py) to your netId. (Do not include your first name, last name ... or any extra character)
# 2) To run the program type the following statement in the command line:  
#       -) python studentNetId.py DNASeq1FilePath DNASeq2FilePath OutputFilePath                                                                   
#    where  DNASeq1FilePath is the path to the file that contains First DNA sequence (e.g. DNASeq1.txt)
#           DNASeq2FilePath is the path to the file that contains Second DNA sequence (e.g. DNASeq2.txt)
#           OutputFilePath is the path that the output is goint to be saved (e.g. out.txt)
# 3) Make sure to replace FirstName, LastName, SectionNumber, NetId in studentInfo with your information
# 4) You may add as many functions as you want to this program
# 5) The core function in your program is DNASeqAlignment function, where it takes three arguments (DNASeq1,DNASeq2,outputPath) and 
#    computes the similarityScore, sequenceAlignment1 and sequenceAlignment2. At the end, the function writes the result to the output file (Do not make any changes to the output section).
# 6) sequenceAlignment1 and sequenceAlignment2 are strings and they are composed of following characters: 'A', 'T', 'G', 'C' and '-', Do not include extra space or any other character in the strings.
# 7) Make sure your program works with at least one of the following python versions: (2.7.9, 2.7.8, 2.6.6)
# 8) Once you have tested your program with one of the versions listed above, assign that version number to pythonVersion in studentInfo function
# 9) Make sure to write enough comments in order to make your code easy to understand. 
# 10) Describe your algorithm in ALGORITHM section below (you may add as many lines as you want).
# 11) To understand the flow of the program consider the following example:
#      0) Let say we have DNASeq1.txt file which contains AACCTGACATCTT and DNASeq2.txt file contains CCAGCGTCAACTT
#      1) If we execute the following command in the command line: -) python studentNetId.py DNASeq1.txt DNASeq2.txt out.txt
#      2) input arguments are parsed       
#      3) studentInfo() function will be executed and the output will be saved in out.txt file
#      4) DNASeqAlignment() function will be called
#      5) At the entry of the DNASeqAlignment function, DNASeq1='AACCTGACATCTT' and DNASeq2='CCAGCGTCAACTT'
#      6) You should compute the sequence alignment of DNASeq1 and DNASeq2. Let say the result is as follows:
#       A A C C T G A C - - - - A T C T T
#       | | | | | | | | | | | | | | | | |
#       - - C C A G - C G T C A A - C T T      
#      7) At the end of the DNASeqAlignment function sequenceAlignment1='AACCTGAC----ATCTT', sequenceAlignment2='--CCAG-CGTCAA-CTT', similarityScore=6.25
#      8) In the output section the result is going to be saved in out.txt file
#-------------------------------------------------------------------------------

# ALGORITHM: 
#
#
#
#
#


import os
import sys
import argparse

def studentInfo():
    pythonVersion = '2.7.9'
    studentFirstName = "Mark"
    studentLastName = "Asfour"
    studentSectionNumber = "2"
    studentNetId = "masfo001"
    info = 'FirstName: ' + studentFirstName + '\n'
    info = info + 'LastName: ' + studentLastName + '\n'
    info = info + 'Section: ' + studentSectionNumber + '\n'
    info = info + 'NetId: ' + studentNetId + '\n'
    info = info + 'Python version: ' + pythonVersion + '\n'
    return info

def biggest(a, b, c):
    Max = a
    if b > Max:
        Max = b
    if c > Max:
        c = Max
    return Max

def diff(x, y, DNASeq1, DNASeq2):
    if (DNASeq1[x] == DNASeq2[y]):
        return 1
    elif ((DNASeq1[x] == 'A' and DNASeq2[y] == 'T') or (DNASeq1[x] == 'T' and DNASeq2[y] == 'A')):
        return -.15
    elif ((DNASeq1[x] == 'G' and DNASeq2[y] == 'C') or (DNASeq1[x] == 'C' and DNASeq2[y] == 'G')):
        return -.15
    else:
        return -.1


def DNASeqAlignment(DNASeq1,DNASeq2,outputPath):
    similarityScore = -1
    sequenceAlignment1 = ''
    sequenceAlignment2 = ''
    #########################################################################################
    # Compute new values for similarityScore and sequenceAlignment1 and
    # sequenceAlignment2 #
    #########################################################################################

    Array = []

    col = len(DNASeq1) + 1              #length of first + 1 for blanks
    row = len(DNASeq2) + 1              #length of second + 1 for blanks

    for x in range(col):                #construct an array
        Array.append([0] * row)
 
    for x in range(col):                #init 0th col
        Array[x][0] = 0 - (.2*x)

    for x in range(row):                #init 0th row
        Array[0][x] = 0 - (.2*x)

    for x in range(1, col):             #fill out the table
        for y in range(1, row):
            remove = Array[x - 1][y] - .2
            insert = Array[x][y - 1] - .2
            match = diff(x - 1, y - 1, DNASeq1, DNASeq2) + Array[x - 1][y - 1]
            Array[x][y] = max(remove, insert, match)
    
    similarityScore = Array[col - 1][row - 1] #score = bottom right of array
    
    a = col - 1
    b = row - 1
    c = []
    while (a != 0 and b != 0):          #find sequence
        if (a != 0):
            i = Array[a][b - 1]         #look above
        if (a != 0 and b != 0):
            j = Array[a - 1][b - 1]     #look diagonal
        if ( b != 0):
            k = Array[a - 1][b]         #look left

        h = max(i, j, k)
        if (h == i):
            b = b - 1
            c.append(1)                 
        elif (h == j):
            a = a - 1
            b = b - 1
            c.append(2)                 
        elif (h == k):
            a = a - 1
            c.append(3)

    seq1 = []
    for x in DNASeq1:
        seq1.append(x)
    seq2 = []
    for x in DNASeq2:
        seq2.append(x)
    for x in range (0, len(c)):         #make subsequence
        if (c[x] == 1):                 #if moved up
            seq2[row - x] = '_'
        elif (c[x] == 2):               #if moved diagonal
            continue 
        elif (c[x] == 3):               #if moved left
            seq1[col - x] = '_'

    sequenceAlignment1 = ''.join(seq1)
    sequenceAlignment2 = ''.join(seq2)

    #for c in range(col):
    #    for r in range(row):
    #        print Array[c][r],
    #    print


    #################################  Output Section  ######################################
    result = "Similarity score: " + str(similarityScore) + '\n'
    result = result + "Sequence alignment1: " + sequenceAlignment1 + '\n'
    result = result + "Sequence alignment2: " + sequenceAlignment2 + '\n'
    writeToFile(outputPath,result)
    
def writeToFile(filePath, content):
    with open(filePath,'a') as file:
        file.writelines(content)

def readFile(filePath):
    logLines = ''
    with open(filePath,'r') as file:
        for logText in file:
            logLines = logLines + logText

    uniqueChars = ''.join(set(logLines))
    for ch in uniqueChars:
        if ch not in {'A','a','C','c','G','g','T','t'}:
            logLines = logLines.replace(ch,'')
    logLines = logLines.upper()
    return logLines

def removeFile(filePath):
    if os.path.isfile(filePath):
        os.remove(filePath)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DNA sequence alignment')
    parser.add_argument('DNASeq1FilePath', type=str, help='Path to the file that contains First DNA sequence')
    parser.add_argument('DNASeq2FilePath', type=str, help='Path to the file that contains Second DNA sequence')
    parser.add_argument('OutputFilePath', type=str, help='Path to the output file')
    args = parser.parse_args()
    DNASeq1 = readFile(args.DNASeq1FilePath)
    DNASeq2 = readFile(args.DNASeq2FilePath)
    outputPath = args.OutputFilePath
    removeFile(outputPath)
    writeToFile(outputPath,studentInfo())
    DNASeqAlignment(DNASeq1,DNASeq2,outputPath)
