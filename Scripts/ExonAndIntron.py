import os
import sys
import getopt
from Bio.Seq import Seq
import re
import copy 


def Usage():
    print "\n Application %s [options] -i <SeqFile> -b >blastfile> \n" \
        "-i     Input Genomic Scaffold File \n" \
        "-b     BlastFile with tab delinieated colulmn in specified format \n" \
        "-o     The output file you are going to write to. DO NOT INCLUDE.ending. All FILES will end with a ,out \n" \
        "'7 qseqid qstart qend sseqid sstart send length pident' | \n" \
        "sort -k 3,3n > BLASTINPUTFILE \n" % (sys.argv[0])

def Main():

    Iflag = None
    Bflag = None
    oflag = None
    Nucflag = None

    try:
        options, other_args = getopt.getopt(sys.argv[1:], "i:b:h:o:", ["help"])

    except getopt.GetoptError:
        print "There was a command parsing error"
        Usage()
        sys.exit(1)

    for option, value in options:
        if option == "-i":
            Iflag = value
        elif option == "-b":
            Bflag = value
        elif option == "-o":
            oflag = value
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
    elif oflag == None:
        print "Need output"
        Usage()
        exit(-1)
    
#Main function calls. Actually run upon correct usage needs.
    ReadThatBlast = BlastFilereader(Bflag)
    ReadThatThing = (GenomeReader(Iflag))
    ExonGetter(ReadThatBlast, ReadThatThing, oflag)
    UniqPairs = UniqueGenePairs(ReadThatBlast)
    GenesAndLocation = counterAuth(UniqPairs, ReadThatBlast)
    LocationReturn(GenesAndLocation,ReadThatThing, oflag)
    Y = IntronLocationRetriever(GenesAndLocation)
    IntronPuller(Y, ReadThatThing, oflag)

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
    # BlastFilereader will return a list of list from a blast file with the
    #above commented 7 header feature from a BLAST table."""
    HitsNStuff = []
    with open(BLastfile, "r") as f:
        for line in f:
            HitsNStuff.append(line.strip("\n").split('\t'))
    return HitsNStuff

def ExonGetter(Blastlist, GenomeFile, OutFileName):
    NewFileName = str(OutFileName) + 'ExonOnly.out'
    Logical = os.path.exists(NewFileName)
    if Logical == True:
        os.remove(NewFileName)
    else: 
        pass 
    with open(NewFileName, 'a') as f: 
        for line in Blastlist:
            GenomicScaffold = line[1]
            gene = line[0]
            FinalString = '>' + gene + GenomicScaffold
            if line[4] < line[5]:
                seq = GenomeFile[GenomicScaffold][int(line[4]) - 1:int(line[5])
                        + 1]
                f.write(FinalString)
                f.write('\n')
                f.write(seq)
                f.write('\n')
            elif line[4] > line[5]:
                Startpost = line[5]
                Negatpost = line[4]
                PulledSeq = GenomeFile[GenomicScaffold][int(Startpost)-1:int(Negatpost)+1]
                NewSeq2 = Seq(PulledSeq).reverse_complement()
                Newseq = str(NewSeq2)
                f.write(FinalString)
                f.write('\n')
                f.write(Newseq)
                f.write('\n')

def ExtronLocationGetter(CutBlastList):
    PosotiveDict = {}
    NegativeDict = {}
    for item in CutBlastList:
        if int(item[4]) < int(item[5]):
            names = (item[0], item[1])
            LocationList = [int(item[4]), int(item[5])]
            PosotiveDict[names] = LocationList
        elif int(item[4]) > int(item[5]):
            names = (item[0], item[1])
            LocationList = [int(item[4]), int(item[5])]
            NegativeDict[names] = LocationList
    return PosotiveDict , NegativeDict

def UniqueGenePairs(BlastFile):
    #The purpose of this function is to take in the blast table and find all
    #unique pairs of Genomics Scaffolds and gene Regions. This section is being
    #refactored in hopes of A. Becoming more readable and B becoming more
    #efficient
    DictofPossible = {}
    GeneGroups = []
    Finalized_gene_Location = []
    for item in BlastFile:
        GenomicScaffold = item[0]
        GeneFound = item[1]
        ComboList = [GenomicScaffold, GeneFound]
        GeneGroups.append(ComboList)
    UniquePairs = set(map(tuple, GeneGroups))
    for item in UniquePairs:
        if item not in DictofPossible:
            DictofPossible[item] = [['Posotive'],['Negative']]
    return DictofPossible

def counterAuth(Uniques, BlastFile):
    for item in BlastFile:
        GenomicScaffold = item[0]
        GeneFound = item[1]
        ComboList = tuple([GenomicScaffold, GeneFound])
        Start = item[4]
        End = item[5]
        if Start < End and ComboList in Uniques:
            Uniques[ComboList][0].append(Start)
            Uniques[ComboList][0].append(End)
        elif Start > End and ComboList in Uniques:
            Uniques[ComboList][1].append(Start)
            Uniques[ComboList][1].append(End)
        else:
            "Something You wrote is fucked"
    return Uniques



def IntronLocationRetriever(LocationInDictFormat):
    ClearedDict = {}
    for item, vale in LocationInDictFormat.iteritems():
        ReDone = []
        for thing in vale:
            if len(thing) >= 4:
                LengthOFSeq = len(thing) - 1
                if LengthOFSeq % 2 == 0:
                    Z = IntronFinder(thing)
                    ClearedDict[item] = Z
                else:
                    pass
    return ClearedDict

def IntronFinder(List1):
    Strand = List1.pop(0)
    Z = len(List1) - 1
    RemoveSeq = set([0,int(Z)])
    Newlist = [v for i, v in enumerate(List1) if i not in RemoveSeq]
    MakeInt = [int(x) for x in Newlist]
    TestSet = set(MakeInt)
    if len(TestSet) != len(MakeInt):
        Y = DuplicateFinder(List1)
        CleanedList = []
        for item in List1:
            if item == Y:
                pass
            else:
                CleanedList.append(item)
        NewLen = len(CleanedList) - 1 
        RemoveSeq2 = set([0,int(NewLen)])
        Newlist2 = [v for i, v in enumerate(CleanedList) if i not in RemoveSeq2]
        RevisedList = [int(x) for x in Newlist2]
        SortedRevised = sorted(RevisedList)
        SortedRevised.insert(0, Strand)
        return SortedRevised
    else:
        SortedNewList = sorted(MakeInt)
        SortedNewList.insert(0, Strand)
        return SortedNewList

def DuplicateFinder(List1):
    uniques = set(List1)
    for item in uniques:
        Count = List1.count(item)
        if Count > 1:
            return item
        else:
            pass

def IntronPuller(DictionaryThing, GenomeFile, OutFileName):
    NewFileName = str(OutFileName) + 'Intron.out'
    Logical = os.path.exists(NewFileName)
    if Logical == True:
        os.remove(NewFileName)
    else:
        pass
    with open(NewFileName, 'a') as f:
        for header, locs in DictionaryThing.iteritems():
            Scaffheader = header[1]
            if locs[0] == 'Posotive':
                Direction = locs.pop(0)
                GenomicSeq = GenomeFile[Scaffheader]
                if len(locs) == 2:
                    IntronicSeq = GenomicSeq[int(locs[0]) - 1:int(locs[1]) + 1]
                    Structure = ('>' + header[0] + " " + header[1] + " " +
                                str(locs[0]) + " " + str(locs[1]) + " " + '+')
                    f.write(Structure)
                    f.write("\n")
                    f.write(str(IntronicSeq))
                    f.write("\n")
                elif len(locs) > 2:
                    n = 2
                    Z = [locs[i:i+n] for i in xrange(0, len(locs), n)]
                    for i in Z:
                        start = int(i[0])
                        end = int(i[1])
                        IntronicSeq = GenomicSeq[start - 1:end + 1]
                        Structure = ('>' + header[0] + " " + header[1] + " " +
                                str(start) + " " + str(end) + " " + '+')
                        f.write(Structure)
                        f.write("\n")
                        f.write(str(IntronicSeq))
                        f.write("\n")
            elif locs[0] == 'Negative':
                Direction = locs.pop(0)
                GenomicSeq = GenomeFile[Scaffheader]
                if len(locs) == 2:
                    IntronicSeq = GenomicSeq[int(locs[0]) - 1:int(locs[1]) + 1]
                    CorrectedSeq = Seq(IntronicSeq).reverse_complement()
                    Structure = ('>' + header[0] + " " + header[1] + " " +
                                str(locs[0]) + " " + str(locs[1]) + " " + '-')
                    f.write(Structure)
                    f.write("\n")
                    f.write(str(CorrectedSeq))
                    f.write("\n")
                elif len(locs) > 2:
                    n = 2
                    Z = [locs[i:i+n] for i in xrange(0, len(locs), n)]
                    for i in Z:
                        start = int(i[0])
                        end = int(i[1])
                        IntronicSeq = GenomicSeq[start - 1:end + 1]
                        CorrectedSeq = Seq(IntronicSeq).reverse_complement()
                        Structure = ('>' + header[0] + " " + header[1] + " " +
                                str(start) + " " + str(end) + " " + '-')
                        f.write(Structure)
                        f.write("\n")
                        f.write(str(CorrectedSeq))
                        f.write("\n")
            else:
                raise Exception("MAJOR ERROR, SEQUENCE HAS NO Polarity ID")

def LocationReturn(DictionaryofNonsense, GenomeFile, OutFileName):
    #This function will take ina  dictionary with a location data in the form
    #of a list with posotive/negative as the first posotion. Will process each
    #of these list and check for same strand, as well as distance being
    #reasonable
    NewFileName2 = str(OutFileName) + 'CombinedExonandIntron.out'
    Logical = os.path.exists(NewFileName2)
    Norp = copy.deepcopy(DictionaryofNonsense)
    if Logical == True:
        os.remove(NewFileName2)
    else: 
        pass 
    with open(NewFileName2, "a") as r:
        for hit, locations in Norp.iteritems():
            ScaffHeader = hit[1]
            for item in locations:
                if len(item) == 3: 
                    pass
                elif len(item) > 3:
                    if item[0] == "Posotive":
                        item.pop(0)
                        item = [ int(x) for x in item ]
                        Max1 = max(item)
                        Min1 = min(item)
                        GenomeSeq = GenomeFile[ScaffHeader]
                        Sequence = GenomeSeq[Min1-1:Max1+1]
                        Structure = ('>' + hit[0] + " " + hit[1] + " " +
                                str(Min1) + " " + str(Max1) + " " + '+')
                        r.write(Structure)
                        r.write("\n")
                        r.write(Sequence)
                        r.write("\n")
                    elif item[0] == "Negative":
                        item.pop(0)
                        item = [ int(x) for x in item ]
                        Max1 = max(item)
                        Min1 = min(item)
                        GenomeSeq = GenomeFile[ScaffHeader]
                        Segment = GenomeSeq[Min1-1:Max1+1]
                        ProperSeq = Seq(Segment).reverse_complement()
                        Structure = ('>' + hit[0] + " " + hit[1] + " " +
                                str(Min1) + " " + str(Max1) + " " + '-')
                        r.write(Structure)
                        r.write("\n")
                        r.write(str(ProperSeq))
                        r.write("\n")


if __name__ == '__main__':
    Main()


