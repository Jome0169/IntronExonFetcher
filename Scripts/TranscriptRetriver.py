import os
import sys
import getopt
from Bio.Seq import Seq
import re
import copy 
"""
This program takes in a trinity assebly and a BLAST file which used the
Trinity.fasta file as referance. This program then goes through, and finds all
uniq isoforms that hit to the BLAST> If multiple thing hit to the SAME isoform,
then the individual with the lowest EVALUE is used. This can be seen in the
header information. 

"""

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

#Main function calls. Actually run upon correct usage needs.

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
    Reads in the BLAST table and looks for uniq isoforms. If the item is
    already within the dicitonary storing names, it appends the blast info to a
    list within the VALUE of the dictionary. Should return a dict with a list
    as the value pair with ONE blast output. (lowest e score)
    
    """ 
    Isoforms = {}
    for item in BlastLocationList: # Adds all Headers 
        if item[1] not in Isoforms:
            Isoforms[str(item[1])] = []
    for item2 in BlastLocationList: # Adds the rest of the sequence
        if item2[1] in Isoforms:
            Isoforms[item2[1]].append(item2)
        else:
            pass
    for key, value in Isoforms.iteritems():
        value.sort(key=lambda x: float(x[7])) # Sorts all values based on Escore 
        del value[1:] #Deletes all but the LOWEST evalue (Closest match)
    return Isoforms


def TrinityReader(GenomeFile):
    """
    Reads in the TRINITY.fasta file into a Dict like object, with the isoforms
    being the key and the Sequence associated with them being the Value. 

    This Function expects the '>' to be in the file to indicate header
    Sequences
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
    Takes in Nonesense, and writes it. Makes sure gene is in Trinity and then
    writes corresponding Sequnece to output File
    """
    TotalStrings = []
    for item in ImportantScaff:
        for thing in Trinity:
            if item in thing:
                if int(ImportantScaff[item][0][4]) <= int(ImportantScaff[item][0][5]):
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
                    FinalizedHeader += '+'
                    StringCreation.append(FinalizedHeader)
                    StringCreation.append(Trinity[thing])
                    TotalStrings.append(StringCreation)
                elif int(ImportantScaff[item][0][4]) >= int(ImportantScaff[item][0][4]):
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
                    FinalizedHeader += '-'
                    StringCreation.append(FinalizedHeader)
                    StringCreation.append(Trinity[thing])
                    TotalStrings.append(StringCreation)
    return TotalStrings


def FileWriter(outputflag, stringarg):
    with open(outputflag, 'a') as f:
        for nested in stringarg:
            for item in nested:
                f.write(item)
                f.write('\n')



if __name__ == '__main__':
    Main()
