import os 
from sys import argv

First, Second, Third= argv

Headers = {}

def GenomeReader(GenomeFile):
    GenomeScaffolds = {}
    with open(GenomeFile, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if line in GenomeScaffolds:
                    NamedSeq = line
                else:
                    NamedSeq = line
                    GenomeScaffolds[NamedSeq] = ""
            else:
                GenomeScaffolds[NamedSeq] += line
        return GenomeScaffolds

def DictWriter(Dictionary):
    with open (Third, 'a') as R:
        for item, thing in Dictionary.iteritems():
            R.write(item)
            R.write('\n')
            R.write(thing)
            R.write('\n')

Z = GenomeReader(Second)
DictWriter(Z)





