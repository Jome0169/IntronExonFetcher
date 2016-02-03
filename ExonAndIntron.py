import os
import sys
import getopt
from Bio.Seq import Seq
import re



def Usage():
    print "\n Application %s [options] -i <SeqFile> -b >blastfile> \n" \
        "-i     Input Genomic Scaffold File \n" \
        "-b     BlastFile with tab delinieated colulmn in specified format \n" \
        "'7 qseqid qstart qend sseqid sstart send length pident' | \n" \
        "sort -k 3,3n > BLASTINPUTFILE \n" % (sys.argv[0])


def Main():

    Iflag = None
    Bflag = None

    try:
        options, other_args = getopt.getopt(sys.argv[1:], "i:b:", ["help"])

    except getopt.GetoptError:
        print "There was a command parsing error"
        Usage()
        sys.exit(1)

    for option, value in options:
        if option == "-i":
            Iflag = value
        elif option == "-b":
            Bflag = value
        elif option == "--help":
            Usage()
        else:
            print "Unhandeled options %s %s" % (options)

    if Iflag == None:
        print "Need a Genomic Scaffold input"
        Usage()
        exit(-1)
    elif Bflag == None:
        print "Need a blast file to referance"
        Usage()
        exit(-1)

        print "Processing Fles %s %s" % (Iflag, Bflag)

    #ReadThatThing = GenomeReader(Iflag)
    #ReadThatBlast = BlastFilereader(Bflag)
    #GetExonAndIntro = ExonAndIntronGetter(ReadThatBlast)
    #DictionaryReaderAndExtractor(ReadThatThing, ReadThatBlast)
    #DictionaryReaderAndExtractor2(ReadThatThing, GetExonAndIntro)
    GeneRipperAndWriter(Bflag)
    


def GenomeReader(GenomeFile):
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



def BlastFilereader(BLastfile):
    """ BlastFilereader will return a list of list from a blast file with the
    above commented 7 header feature from a BLAST table."""
    HitsNStuff = []
    with open(BLastfile, "r") as f:
        for line in f:
            HitsNStuff.append(line.split('\t'))
        return HitsNStuff
    
    
def ExonAndIntronGetter(FilesHitsFromBlast):
    FilthyList = []
    Iterlist = xrange(0,len(FilesHitsFromBlast))
    counter = 1
    trillList = []
    for item in FilesHitsFromBlast:
            for number in Iterlist:
                if (item[0] == FilesHitsFromBlast[number][0] and item[3] ==
                FilesHitsFromBlast
                [number][3] and item[1] != FilesHitsFromBlast[number][1]):
                    trillList.append([item[0], item[3]])
                else:
                    pass
    UniqGeneGroups = set(map(tuple, trillList))
    for item in UniqGeneGroups:
        WEGONNAGETIT = []
        for line in FilesHitsFromBlast:
            if item[0] == line[0] and item[1] == line[3]:
                WEGONNAGETIT.append(line)
        NumList = []
        NumList2 = []
        for item in WEGONNAGETIT:
            NumList.append(int(item[1]))
            NumList2.append(int(item[2]))
        Start = min(NumList)
        End = max(NumList2)
        Thing = WEGONNAGETIT[0]
        Thing[1] = int(Start)
        Thing[2] = int(End)
        FilthyList.append(Thing)
    return FilthyList


def DictionaryReaderAndExtractor(DictionaryForUse, BlastFileOutput):
    for item in BlastFileOutput:
        if item[0] in DictionaryForUse:
            if item[4] < item[5]:
                SequenceMix = ""
                StartPos = int(item[1]) - 1
                EndPos = int(item[2])
                SequenceMix += DictionaryForUse.get(item[0])
                NewSeq = SequenceMix[StartPos:EndPos]
                FastaHeaderName = re.search("(scaffold[1-9]{1,8})\w",str(item[0]))
                ScaffoldHeader = FastaHeaderName.group(0)
                NewheaderName = ">" + str(ScaffoldHeader) + "_" + str(item[3]) 
                ExonFileWriter(NewheaderName)
                ExonFileWriter('\n')
                ExonFileWriter(str(NewSeq))
                ExonFileWriter('\n')
            elif item[5] < item[4]:
                SequenceMix = ""
                StartPos = int(item[1]) - 1
                EndPos = int(item[2])
                SequenceMix += DictionaryForUse.get(item[0])
                PulledSeq = SequenceMix[StartPos:EndPos]
                NewSeq = Seq(PulledSeq).reverse_complement()
                FastaHeaderName = re.search("(scaffold[1-9]{1,8})\w",str(item[0]))
                ScaffoldHeader = FastaHeaderName.group(0)
                NewheaderName = ">" + str(ScaffoldHeader) + "_" + str(item[3]) 
                ExonFileWriter(NewheaderName)
                ExonFileWriter('\n')
                ExonFileWriter(str(NewSeq))
                ExonFileWriter('\n')

def DictionaryReaderAndExtractor2(DictionaryForUse, BlastFileOutput):
    for item in BlastFileOutput:
        if item[0] in DictionaryForUse:
            if item[4] < item[5]:
                SequenceMix = ""
                StartPos = int(item[1]) - 1
                EndPos = int(item[2])
                SequenceMix += DictionaryForUse.get(item[0])
                NewSeq = SequenceMix[StartPos:EndPos]
                FastaHeaderName = re.search("(scaffold[1-9]{1,8})\w",str(item[0]))
                ScaffoldHeader = FastaHeaderName.group(0)
                NewheaderName = ">" + str(ScaffoldHeader) + "_" + str(item[3]) 
                IntronExonFileWriter(NewheaderName)
                IntronExonFileWriter('\n')
                IntronExonFileWriter(str(NewSeq))
                IntronExonFileWriter('\n')
            elif item[5] < item[4]:
                SequenceMix = ""
                StartPos = int(item[1]) - 1
                EndPos = int(item[2])
                SequenceMix += DictionaryForUse.get(item[0])
                PulledSeq = SequenceMix[StartPos:EndPos]
                NewSeq = Seq(PulledSeq).reverse_complement()
                FastaHeaderName = re.search("(scaffold[1-9]{1,8})\w",str(item[0]))
                ScaffoldHeader = FastaHeaderName.group(0)
                NewheaderName = ">" + str(ScaffoldHeader) + "_" + str(item[3]) 
                IntronExonFileWriter(NewheaderName)
                IntronExonFileWriter('\n')
                IntronExonFileWriter(str(NewSeq))
                IntronExonFileWriter('\n')


def GeneRipperAndWriter(BlastfiletoRead):
    AllGenesPresent = []
    with open(BlastfiletoRead, 'r') as f:
        for item in f:
            x = item.split("\t")
            AllGenesPresent.append(x[3])
        OnlyUniqGenes = set(AllGenesPresent)
        ListReConvert = list(OnlyUniqGenes)
        print ListReConvert


def ExonFileWriter(ArgumentToWrite):
    with open("ExonFile.fasta", 'a') as f:
        f.write(ArgumentToWrite)

def IntronExonFileWriter(ArgumentToWrite):
    with open("IntronAndExonFile.fasta", 'a') as f:
        f.write(ArgumentToWrite)

if __name__ == '__main__':
    Main()








