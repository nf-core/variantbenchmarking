#!/usr/bin/env python
import re
import sys
import os

#make mates dictionary given list input of non-comment lines
def makeMateDict(m):
    d = {}
    for index1, line1 in enumerate(m):
        id1 = line1.split('\t')[2]
        numMate = re.search(r':(\d)',id1).group(1)
        origId = re.search(r'(\d+):',id1).group(1)
        if int(numMate) == 1:
            for index2, line2 in enumerate(m):
                #never start from beginning of file
                if index2 <= index1:
                    continue
                # print str(index1) + " : " + str(index2)
                id2 = line2.split('\t')[2]
                duplicateId = re.search(r'(\d+):',id2).group(1)
                duplicateNumMate = re.search(r':(\d)',id2).group(1)
                if duplicateId == origId and int(duplicateNumMate) == 2:
                    d[line1] = line2
                    break
    return d

def classify(line, ALT_INDEX, mdict):
    #get alt, chrom1, chrom2, position (pos), id, old SVTYPE (should be BND if virgin svaba vcf) from line
    s = line.split("\t")
    alt = s[ALT_INDEX]
    chrom1 = s[0]
    pos = int(s[1])
    id=s[2]

    if int(re.search(r':(\d)',id).group(1)) != 1:
        return "NONE"

    mateLine = mdict[line].split('\t')
    mateChrom = mateLine[0]
    mateAlt = mateLine[ALT_INDEX]

    oldType = re.search(r'SVTYPE=(.+?)(\s+?|:)',line).group(1)

    # get new type
    if oldType == 'BND' and chrom1 == mateChrom:
        INV_PATTERN_1 = re.compile(r'\D\].+\]')
        INV_PATTERN_2 = re.compile(r'\[.+\[\D')
        if INV_PATTERN_1.match(alt) and INV_PATTERN_1.match(mateAlt):
            return "INV"
        if INV_PATTERN_2.match(alt) and INV_PATTERN_2.match(mateAlt):
            return "INV"

        # DEL
        DEL_PATTERN_THIS = re.compile(r'\D\[.+\[')
        DEL_PATTERN_MATE = re.compile(r'\].+\]\D')
        if DEL_PATTERN_THIS.match(alt) and DEL_PATTERN_MATE.match(mateAlt):
            return "DEL"

        # INS
        INS_PATTERN_THIS = re.compile(r'\D\].+\]')
        INS_PATTERN_MATE = re.compile(r'\[.+\[\D')
        if INS_PATTERN_THIS.match(alt) and INS_PATTERN_MATE.match(mateAlt):
            return "DUP/INS"

    return 'BND'

if __name__ == "__main__":
    file = sys.argv[1]
    if not os.path.exists(file):
        raise IOError(file)
    alt_index = -1
    #generate mate:mate dictionary
    #load file into ram
    vcf_file=[]
    with open (file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            vcf_file.append(line)
    matesDict = makeMateDict(vcf_file)
    with open(file, "r") as f:
        for line in f:
            # print comments
            if line.startswith("##"):
                sys.stdout.write(line)
                continue
            # header contains indexes
            if line.startswith('#'):
                split = line.split("\t")
                for index, val in enumerate(split):
                    if val == "ALT":
                        alt_index = index
                        break
                sys.stdout.write(line)
                continue
            if alt_index == -1:
                print "ERROR: NO ALT INDEX FOUND"
                exit(1)
            newType = classify(line, alt_index, matesDict)
            if newType != "NONE":
                newLine = re.sub(r'SVTYPE=BND',"SVTYPE="+newType,line)
                sys.stdout.write(newLine)