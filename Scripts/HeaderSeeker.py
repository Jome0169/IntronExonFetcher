import os
from sys import argv

script, file1 = argv

with open(file1, 'r') as f:
    for line in f:
        if line.startswith('@') and ' ' in line:
            print line 
