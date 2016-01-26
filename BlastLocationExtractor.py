import os
import sys
from Bio.Seq import Seq



HitsNStuff = []
FilthyList = []



def BlastFilereader(BLastfile):
    with open(BLastfile, "r") as f:
        for line in f:
            HitsNStuff.append(line.split('\t'))
        Iterlist = xrange(0,len(HitsNStuff))
        counter = 1
        trillList = []
        for item in HitsNStuff:
            for number in Iterlist:
                if (item[0] == HitsNStuff[number][0] and item[3] ==
                HitsNStuff[number][3] and item[1] != HitsNStuff[number][1]):
                    trillList.append([item[0], item[3]])
                else:
                    pass
        UniqGeneGroups = set(map(tuple, trillList))
        for item in UniqGeneGroups:
            WEGONNAGETIT = []
            for line in HitsNStuff:
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

GenomeScaffolds = {}

def GenomeReader(GenomeFile):
    key = []
    value = []
    with open(GenomeFile, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                NamedSeq = line.replace('>', '')
                key.append(NamedSeq)
                GenomeScaffolds[NamedSeq] = ""
            else:
                GenomeScaffolds[NamedSeq] += line


def DictionaryReaderAndExtractor(DictionaryForUse, BlastFileOutput):
    SequenceForFile = []
    for item in BlastFileOutput:
        if item[0] in DictionaryForUse:
            if item[4] < item[5]:
                SequenceMix = ""
                StartPos = int(item[1]) - 1
                EndPos = int(item[2])
                SequenceMix += DictionaryForUse.get(item[0])
                NewSeq = SequenceMix[StartPos:EndPos]
                ExonFileWriter(str(item[0]))
                ExonFileWriter(str(item[3]))
                ExonFileWriter("Posotive Strand")
                ExonFileWriter(str(NewSeq))
                ExonFileWriter('\n')
            elif item[5] < item[4]:
                SequenceMix = ""
                StartPos = int(item[1]) - 1
                EndPos = int(item[2])
                SequenceMix += DictionaryForUse.get(item[0])
                PulledSeq = SequenceMix[StartPos:EndPos]
                NewSeq = Seq(PulledSeq).reverse_complement()
                ExonFileWriter(str(item[0]))
                ExonFileWriter(str(item[3]))
                ExonFileWriter("Negative Strand")
                ExonFileWriter(str(NewSeq))
                ExonFileWriter('\n')



def DictionaryReaderAndExtractor2(DictionaryForUse, BlastFileOutput):
    SequenceForFile = []
    for item in BlastFileOutput:
        if item[0] in DictionaryForUse:
            if item[4] < item[5]:
                SequenceMix = ""
                StartPos = int(item[1]) - 1
                EndPos = int(item[2])
                SequenceMix += DictionaryForUse.get(item[0])
                NewSeq = SequenceMix[StartPos:EndPos]
                IntronExonFileWriter(str(item[0]))
                IntronExonFileWriter(str(item[3]))
                IntronExonFileWriter("Posotive Strand")
                IntronExonFileWriter(str(NewSeq))
                IntronExonFileWriter('\n')
            elif item[5] < item[4]:
                SequenceMix = ""
                StartPos = int(item[1]) - 1
                EndPos = int(item[2])
                SequenceMix += DictionaryForUse.get(item[0])
                PulledSeq = SequenceMix[StartPos:EndPos]
                NewSeq = Seq(PulledSeq).reverse_complement()
                IntronExonFileWriter(str(item[0]))
                IntronExonFileWriter(str(item[3]))
                IntronExonFileWriter("Negative Strand")
                IntronExonFileWriter(str(NewSeq))
                IntronExonFileWriter('\n')


 


def ExonFileWriter(ArgumentToWrite):
    with open("ExonFile.pabs", 'a') as f:
        f.write(ArgumentToWrite)
        f.write('\n')

def IntronExonFileWriter(ArgumentToWrite):
    with open("IntronAndExonFile.pabs", 'a') as f:
        f.write(ArgumentToWrite)
        f.write('\n')


BlastFilereader(sys.argv[1])
GenomeReader(sys.argv[2])
DictionaryReaderAndExtractor(GenomeScaffolds, HitsNStuff)
DictionaryReaderAndExtractor2(GenomeScaffolds, FilthyList)
