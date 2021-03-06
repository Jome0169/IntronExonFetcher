import os, sys
import getopt
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

global Start
global End

Start = ['ATG']
End = ['TAA', 'TAG', 'TGA']

def GenomeReader(GenomeFile):
    """Reads in a Multifasta file and returns a dictionary. The dictionary works
    in such a way as the key is equal to the header and the value is the
    Seqeunece"""
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
        return GenomeScaffolds  # Returns a Dictionary object

class SequencePossible(object):
    """
    The purpose of this class is to take in a sequence object mush like
    biopython and return a host of different reading frames that could be the
    protein coding region. 

    It should be noted that this class does not actually check for strand
    polarity. Instead that must nbe done before this class is called on the
    sequecne object. 
    """
    
    def __init__(self, Seqq):
        """ 
        Intialize two list to store protein seq as well as as the actual
        sequence itself and the sequence length
        """
        self.Seqq = Seqq
        self.length = len(Seqq)
        self.FinalMaxProt = []
        self.FinalMacNuc = []

    def RevComp(self):
        """
        Uses biopython to take reverse complement of nucleotoide sequence if it
        was found to be negative in previous steps
        """
        Y = Seq(self.Seqq).reverse_complement()
        return self.Seqq

    def SlidingWindow(self):
        """
        This function looks up every possible start and stop sequence within
        the given Sequence by looking for ATG and any of the three possible
        stop sequences. Using a combination of all of them to add all protein
        sequences to the above self.Proteins function. Creates a list of list
        of possible proteins 
        """
        ProteinStorage = []
        NucStorage = []
        SeqToEdit = self.Seqq
        try:
            iter(SeqToEdit)
        except:
            raise Exception("Must Be iterable Seq")
        Z = xrange(0,self.length)
        Starts = []
        Ends = []
        for i in Z:
            if SeqToEdit[i:i+3] in Start:
                #print i, i+3
                #print SeqToEdit[i:i+3]
                Newrane = xrange(i,self.length, 3)
                for item in Newrane:
                    Stogare = str(SeqToEdit[item:item+3])
                    if Stogare in End:
                        DNA = SeqToEdit[i:item+3]
                        NucStorage.append(DNA)
                        MessengerRNA = Seq(DNA).transcribe()
                        Protein = MessengerRNA.translate()
                        ProteinStorage.append(Protein)
                        break
        if not ProteinStorage:
            pass
        else:
            LongestProtein = max(ProteinStorage, key=len)
            LongestNuc = max(NucStorage, key=len)
            self.FinalMaxProt.append(str(LongestProtein))
            self.FinalMacNuc.append(str(LongestNuc))


    def FileWriter(self, FileToWrite, Headername):
        ProtToWrite = self.FinalMaxProt
        NucToWrite = self.FinalMacNuc
        with open(FileToWrite, 'a') as f:
            StringHeader = '>' + str(Headername) 
            for item in ProtToWrite:
                if len(item) > 110:
                    X = StringHeader + str(len(item))
                    f.write(X)
                    f.write('\n')
                    f.write(str(item))
                    f.write('\n')
                else:
                    pass
        Nuc = str(FileToWrite) + 'NucSeq'
        with open(Nuc, 'a') as f:
            StringHeader = '>' + str(Headername)
            for item in NucToWrite:
                f.write(StringHeader)
                f.write('\n')
                f.write(str(item))
                f.write('\n')



def Usage():
    print "\n Application %s [options] -i <SeqFile> -o <Output>\n" \
        " NOTE: The INPUT FILE must be of the same format as a file which has \n" \
        "been outputted from NewExonIntron.py script. Can be found on my githun\n" \
        "-i     Input Genomic File \n" \
        "-o     Thr Output file. Wil be a translated Protein file \n"  % (sys.argv[0])

def Main():

    Iflag = None
    output = None

    try:
        options, other_args = getopt.getopt(sys.argv[1:], "i:o:h:", ["help"])

    except getopt.GetoptError:
        print "There was a command parsing error"
        Usage()
        sys.exit(1)

    for option, value in options:
        if option == "-i":
            Iflag = value
        elif option == "-o":
            output = value
        elif option == "--help":
            Usage()
        elif option == '-h':
            Usage()
        else:
            print "Unhandeled options %s %s" % (options)

    if Iflag == None:
        print "Need a Genomic Scaffold input"
        Usage()
        exit(-1)
    elif output == None:
        print "Need an Output File to write to "
        Usage()
        exit(-1)
    
    try:
        X = GenomeReader(Iflag)
        for item, value in X.iteritems():
            if item.endswith('+'):
                Z = SequencePossible(value)
                Z.SlidingWindow()
                Z.FileWriter(output, item)
            elif item.endswith('-'):
                Z = SequencePossible(value)
                Z.RevComp()
                Z.SlidingWindow()
                Z.FileWriter(output, item)

    except IOError:
        print "File was not located in this location "


if __name__ == '__main__':
    Main()








