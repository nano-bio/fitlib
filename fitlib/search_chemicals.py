from fitlib.chemlib import ChemicalObject
import fitlib.helplib as hl
import sys

import argparse

#this programme reads a list of chemicals and tries to find InChI and Chemspider ID for each entry.
#be aware that you need a chemspider_token.txt in the directory for the app to work
#the chemspider_token.txt should only contain the token (available online for free)

#using the argparse module to make use of command line options
parser = argparse.ArgumentParser(description="Read a list of chemicals and output the same list completed with InChI and Chemspider IDs (tab separated)")

parser.add_argument('--version', action='version', version='r1')
parser.add_argument('filename', help='Specify a filename to read.')

#parse it
args = parser.parse_args()


try:
    liste = hl.openfile(args.filename)
except IOError:
    sys.exit('Could not read file')
    
try:
    outputliste = open(args.filename + '_output.txt', 'w')
except:
    sys.exit('Could not open output file')

for line in liste:
    line = line.rstrip('\r\n')
    print('Searching for: ' + line)
    try:
        compound = ChemicalObject(name = line)
    except IOError:
        sys.exit('No chemspider token given. Aborting.')
        
    compound.complete()
    outputliste.write(line.decode('utf8', 'replace').encode('ascii', 'replace') + '\t' + compound.inchi + '\t' + compound.inchikey + '\t' + str(compound.csid) + '\n')

liste.close()
outputliste.close()
