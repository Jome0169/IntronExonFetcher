import os
import sys
import getopt
from Bio.Seq import Seq
import re
import copy 


def Usage():
    print "\n Application %s [options] -i <SeqFile> -b >blastfile> \n" \
        "-b     BlastFile with tab delinieated colulmn in specified format \n" \
        "-t     A trinity.fasta file from a trinity transcriptome run \n" \
        "-o     The output file you are going to write to. DO NOT INCLUDE.ending. All FILES will end with a ,out \n" \
        "'7 qseqid qstart qend sseqid sstart send length pident' | \n" \
        "sort -k 3,3n > BLASTINPUTFILE \n" % (sys.argv[0])

def Main():

    Bflag = None
    TrinFlag = None
    oflag = None

    try:
        options, other_args = getopt.getopt(sys.argv[1:], "t:b:h:o:", ["help"])

    except getopt.GetoptError:
        print "There was a command parsing error"
        Usage()
        sys.exit(1)

    for option, value in options:
        if option == "-b":
            Bflag = value
        elif option == "-o":
            oflag = value
        elif option == '-t':
            TrinFlag = value
        elif option == "--help":
            Usage()
        else:
            print "Unhandeled options %s %s" % (options)

    if Bflag == None:
        print "Need a blast file to referance"
        Usage()
        exit(-1)
    elif oflag == None:
        print "Need output"
        Usage()
        exit(-1)
    elif TrinFlag == None:
        print "Need Trinity Input"
        Usage()
        exit(-1)

#Main function calls. Actually run upon correct usage needs. All functions here
#usually return a specified Data type

    BlastLocation = ReadBlastFile(Bflag)
    ImportantIsos = IsoFormDict(BlastLocation)
    TrinityFile = TrinityReader(TrinFlag)
    IsoformsWithSeq = IsoFormRetriever(TrinityFile, ImportantIsos)
    FileWriter(oflag, IsoformsWithSeq)




def ReadBlastFile(BlastFile):
    """ So it appears that withi each of these BLAST tables, a few isoforms are
    creating a majority of the hits here, So what we need to do, and what this
    function does is it takes in the input of the blast file, finds all UNIQ
    isoforrms, and then finds the blast line with the LOWEST e score from this
    set and assumes that is the CORRECT ONE.
    """
    AllLines = []
    with open(BlastFile, 'r') as f:
        for line in f:
            CleanedLine = line.strip('\n').split('\t')
            AllLines.append(CleanedLine)
    return AllLines


def IsoFormDict(BlastLocationList):
    """
    Creates a dictionary of all BLAST hits that matches a certain isoform. Then
    sorts these based off the percent ID and then takes the lowest hitting
    escore (most similar to sequence.) BAsed off some of the preliminary
    outputs however, it is probably going to be best to sort out low sequence
    identities in the beggining. 

    """
    Isoforms = {}
    for item in BlastLocationList:
        if item[1] not in Isoforms:
            Isoforms[str(item[1])] = []
    for item2 in BlastLocationList:
        if item2[1] in Isoforms:
            Isoforms[item2[1]].append(item2)
        else:
            pass
    for key, value in Isoforms.iteritems():
        value.sort(key=lambda x: float(x[7]))
        del value[1:]
    return Isoforms


def TrinityReader(GenomeFile):
    """
    Reads in the trinity file and creates a dictionary using the header and
    isoform information as the key. The Nucleotide Sequence is the value. 
    """

    GenomeScaffolds = {}
    key = []
    with open(GenomeFile, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                NamedSeq = line.replace('>', '')
                key.append(NamedSeq)
                GenomeScaffolds[NamedSeq] = ""
            else:
                GenomeScaffolds[NamedSeq] += line
        return GenomeScaffolds

def IsoFormRetriever(Trinity, ImportantScaff):
    """
    Matches isoforms from the BLAST table with isoforms from the trinity
    assembly. When BLAST isoforms matches it creates a header and printing the
    ISOFORM, TERPENPOIDMATCH, LENGTH, and PIDENT
    """

    TotalStrings = []
    for item in ImportantScaff:
        for thing in Trinity:
            if item in thing:
                StringCreation = []
                thing.split(' ')[1]
                FinalizedHeader = '>'
                FinalizedHeader += str(item)
                FinalizedHeader += '__'
                FinalizedHeader += str(ImportantScaff[item][0][0])
                FinalizedHeader += '__'
                FinalizedHeader += str(thing.split(' ')[1])
                FinalizedHeader += '__'
                FinalizedHeader += str(ImportantScaff[item][0][7])
                StringCreation.append(FinalizedHeader)
                StringCreation.append(Trinity[thing])
                TotalStrings.append(StringCreation)
    return TotalStrings


def FileWriter(outputflag, stringarg):
    """
    Writes a List of List to a file using the outputflag from get opt
    """
    with open(outputflag, 'a') as f:
        for nested in stringarg:
            for item in nested:
                f.write(item)
                f.write('\n')



if __name__ == '__main__':
    Main()
