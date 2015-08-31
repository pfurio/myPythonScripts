#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author: Pedro Furió Tarí

import getopt
import sys
import os.path

class Target:
    def __init__(self, ensID, geneSym, transID, pct):
        self.ensID   = ensID
        self.geneSym = geneSym
        self.transID = transID
        self.pct     = float(pct) if pct != "NULL" else 0.0

    def getPCT(self):
        return self.pct

    def getEnsID(self):
        return self.ensID

    def getGeneSym(self):
        return self.geneSym

    def getTransID(self):
        return self.transID

    def getMiRNA(self):
        return self.miRNA


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hm:t:o:s:", ["help", "mFamily=", "targets=","output=", "speciesID="])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    mfamily = None
    targetF = None
    outputF = None
    speciesid = "9606"

    for o, a in opts:
        if o in ("-h","--help"):
            usage()
            sys.exit()
        elif o in ("-m", "--mFamily"):
            if os.path.isfile(a):
                mfamily = a
            else:
                sys.stderr.write("\nERROR: miR_Family_Info file not recognized.\n")
                usage()
                sys.exit()
        elif o in ("-t", "--targets"):
            if os.path.isfile(a):
                targetF = a
            else:
                sys.stderr.write("\nERROR: Targets file not recognized.\n")
                usage()
                sys.exit()
        elif o in ("-s", "--speciesID"):
            speciesid = a
        elif o in ("-o", "--output"):
            outputF = a
        else:
            assert False, "Unhandled option"

    if mfamily is not None and targetF is not None and outputF is not None:
        run(mfamily, targetF, outputF, speciesid)
    else:
        usage()


def usage():
    print "\nThis script will analyze the two files downloaded from targetScan and will output one line per (mature miRNA - target mRNA)"
    print "\nUsage: python targetScanExtract.py [options] <mandatory>"
    print "Options:"
    print "\t-h, --help:\n\t\t Shows this help"
    print "\t-s, --speciesID\n\t\t Species ID from which we want to extract the pairs (Default: 9606) - http://www.targetscan.org/mmu_61/docs/species.html"
    print "Mandatory:"
    print "\t-m, --mFamily:\n\t\t miR_Family_Info.txt downloaded from targetScan database"
    print "\t-t, --targets:\n\t\t Conserved_Family_Conserved_Targets_Info.txt or Nonconserved_Family_Info.txt from targetScan database"
    print "\t-o, --output:\n\t\t Output file"
    print "\n26/08/2015. Pedro Furió Tarí.\n"


# This function will correct the family IDs given by targetScan when more than one family are found to have 
# identical features
def correctIDs(familiesID):
    myout = []
    lastID = ""
    for family in familiesID:
        f_split = family.split("-")[0]
        if f_split in ["miR", "let"]:
            lastID = f_split
            myout.append(family)
        else:
            myout.append(lastID + "-" + family)
    return myout


def run(mfamily, targetF, outputF, speciesid):

    # 1. We store in a dictionary all the mature miRNAs associated to each miRNA family
    print "Analysing miRNA families file..."
    miR_family = {}

    family_file = open(mfamily, 'r')
    family_file.next() # We skip the first line containing the header
    for line in family_file:
        line_s = line.split()
        
        familiesID = correctIDs(line_s[0].split("/"))
        speciesID  = line_s[2]
        mat_mirna  = line_s[3]

        if speciesID == speciesid:
            for familyID in familiesID:
                if familyID not in miR_family:
                    miR_family[familyID] = []
                miR_family[familyID].append(mat_mirna)
    family_file.close()

    # 2. We check the targets in the target file 
    print "Analysing targets file..."
    target_file = open(targetF, 'r')
    target_file.next() # Skip the first line containing the header
    mytargets = {} # miRNA_family -> [Target(), Target()]
    for line in target_file:
        line_s = line.split()

        familiesID = correctIDs(line_s[0].split("/"))
        ensID      = line_s[1].split(".")[0]
        geneSym    = line_s[2]
        transID    = line_s[3].split(".")[0]
        speciesID  = line_s[4]
        pct        = line_s[10]

        if speciesID == speciesid:
            for familyID in familiesID:
                if familyID in miR_family:
                    if familyID not in mytargets:
                        mytargets[familyID] = []
                    found = False
                    for mytarget in mytargets[familyID]:
                        if transID == mytarget.getTransID():
                            found = True
                            break
                    if found is False:
                        mytargets[familyID].append(Target(ensID, geneSym, transID, pct))

    target_file.close()

    # Get the top 20 targets for each miRNA
    print "Reporting results..."
    myoutput    = open(outputF, 'w')
    myoutput.write ("miRNA_family\tmature_miRNA\tensID\tGeneSymbol\ttransID\tPCT\n")

    for familyID in mytargets:
        sorted_targets = sorted(mytargets[familyID], key = lambda tup:tup.getPCT(), reverse=True)
        # If the pct value is 0, we won't be able to get the top 20, so I will get twice the number of possible targets (although random)
        limit = 20 if sorted_targets[0].getPCT() > 0 else 40
        for mytarget in sorted_targets[:limit]:
            for mat_mirna in miR_family[familyID]:
                myoutput.write(familyID + "\t" + mat_mirna + "\t" + mytarget.getEnsID() + "\t" + mytarget.getGeneSym() + "\t" + 
                    mytarget.getTransID() + "\t" + str(mytarget.getPCT()) + "\n")

    myoutput.close()

if __name__ == "__main__":
    main()

