import sys
import re
import numpy as np
import random
import ete3

import time
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "parsonsc@mit.edu"

# Takes a dictionary of sequences of the form {id1:sequence1, id2:sequence2, ..., idn:sequencen}
# and writes it to the given filename as a fasta file
def writeAlignmentAsFasta(seqDict, outFileName):
    out = open(outFileName, 'w')
    for seq in seqDict:
        out.write('>' + seq + '\n' + seqDict[seq] + '\n')
    out.close()


# Takes in two aligned sequences and returns the proportion of identical gapless sites among the 
# two sequences. If the two sequences have no comparable sites, returns -1.
def percID(seq1, seq2):
    tot = 0
    match = 0
    for i in range(len(seq1)):
        chars = (seq1[i], seq2[i])
        if '-' not in chars:
            tot += 1
            if chars[0] == chars[1]:
                match += 1

    if tot == 0:
        return -1
    return match / tot

# Takes in two aligned sequences and returns False if they are not identical in all gapless sites
# and True if they are.
def checkIdent(seq1, seq2):
    for i in range(len(seq1)):
        if seq1[i] != '-':
            if seq2[i] != '-':
                if seq2[i] != seq1[i]:
                    return False

    return True

def checkIdentGapless(seq1, seq2, sites):
    for i in sites:
        if seq2[i] != '-':
            if seq2[i] != seq1[i]:
                return False

    return True

def checkIdentStoch(seq1, seq2):
    al = list(range(len(seq1)))
    while al:
        i = random.choice(al)
        if seq1[i] != '-':
            if seq2[i] != '-':
                if seq1[i] != seq2[i]:
                    return False
        al.remove(i)

    return True

def checkDuplicateNames(alignmentFileName):
    align = open(alignmentFileName)
    alines = align.readlines()
    align.close()
    names = []
    dups = []
    for line in alines:
        l = line.strip()
        if l[0] == '>':
            if l[1:] in names:
                dups.append(l[1:])
            else:
                names.append(l[1:])

    dups = set(dups)

    return dups



def writeAlignmentAsPhylip(seqDict, outFileName):
    #to be implemented
    pass

# Takes a file of sequences to be removed and an alignment, then returns an alignment
# dictionary with the selected sequences removed. Assumes sequence names in the alignment
# are of the form "culledName_" + anything, but also works if they are simply "culledName"
def cutSeqs(alignmentName, cullFileName):
    cullingFile = open(cullFileName)
    toCullRaw = cullingFile.readlines()
    cullingFile.close()

    toCull = []

    for line in toCullRaw:
        l = line.strip()
        parts = l.split('_')
        if parts[0] in ["WP", "NP", "YP", "XP"]:
            toCull.append('_'.join(parts[0:2]))
        else:
            toCull.append(parts[0])

    seqDict	= fastaToDictionary(alignmentName)

    for bad in toCull:
        if bad in seqDict:
            seqDict.pop(bad)
        else:
            print(bad, "not in alignment??")

    return seqDict

def extractSeqsByTax(alignmentName, keyWord, outFileName, otherExclusions = None):
    NCBItoTax = {}

    IDlist = open("/ncbiToTaxonomy.txt")
    IDlines = IDlist.readlines()
    for line in IDlines:
        l = line.strip()
        parts = l.split()
        #print(parts)
        if len(parts) > 1:
            NCBItoTax[parts[0]] = '_'.join(parts[1:])
            #print('b')
    IDlist.close()

    seqDict = fastaToDictionary(alignmentName)
    
    toPop = []

    if otherExclusions:
        for exc in otherExclusions:
            if exc in seqDict:
                seqDict.pop(exc)
            else:
                print(seq, "not found in alignment.")

    for seq in seqDict:
        if seq in NCBItoTax:
            if keyWord not in NCBItoTax[seq]:
                toPop.append(seq)
        else:
            toPop.append(seq)
            #else:
            #    print(NCBItoTax[seq])

    for seq in toPop:
        seqDict.pop(seq)

    writeAlignmentAsFasta(seqDict, outFileName)

    print(len(toPop), "sequences removed")

def generalizeForReconciliation(treeFileName, taxList, outFileName):
    file = open(treeFileName)
    tree = file.readlines()[0].strip()
    file.close()

    taxCounts = {}
    for tax in taxList:
        taxCounts[tax] = 0

    leaves = getLeafLabels(tree)
    for leaf in leaves:
        simp = ''
        for tax in taxList:
            if tax in leaf:
                if simp == '':
                    simp = tax
                else:
                    print("More than one hit for", leaf + '!')
        if simp == '':
            for i in range(len(taxList)):
                print(str(i) + ". " + taxList[i] + '\n')
            user = input("No good taxonomy found for: " + leaf + '\n')
            simp = taxList[int(user)]
        taxCounts[simp] += 1
        lab = '_' + str(taxCounts[simp])
        tree = tree.replace(leaf, simp + lab)
    #print(taxCounts)
    out = open(outFileName, 'w')
    out.write(tree)
    out.close()



def extractSeqs(alignmentName, pullFileName):
    pullingFile = open(pullFileName)
    toPullRaw = pullingFile.readlines()
    pullingFile.close()

    toPull = []

    for line in toPullRaw:
        l = line.strip()
        toPull.append(l.split('_')[0])

    seqDict = fastaToDictionary(alignmentName)

    for seq in seqDict.copy():
        if seq not in toPull:
            seqDict.pop(seq)
        else:
            toPull.pop(seq)

    for ID in toPull:
        print(ID, "not in alignment??\n")

    return seqDict

# Reads a fasta file given by the provided filename and returns a dictionary of the form
# {id1:sequence1, id2:sequence2, ..., idn:sequencen}
def fastaToDictionary(alignmentName):
    alignFile = open(alignmentName)
    alines = alignFile.readlines()
    alignFile.close

    seqDict = {}
    cur = ""

    for line in alines:
        l = line.strip()
        if l[0] == '>':
            cur = l[1:]
            seqDict[cur] = ""
        else:
            seqDict[cur] += l

    return seqDict

# Returns a list of the leaf labels of a given tree in Newick form. Works with IQTree output.
def getLeafLabels(tree):
    namesTemp = tree.replace("(", "").split(",")
    names = []

    for name in namesTemp:
        parts = name.split(":")
        #names.append(parts[0].replace('>', ''))
        names.append(parts[0])
    return(names)


# Reads a tree of a given name whose leaf labels correspond to NCBI IDs contained within my file
# "/ncbiToTaxonomy.txt". If they are not contained within the file, they are looked up in the
# NCBI database and their taxonomy is added to the file if a match is found. The taxonomy is
# then appended to the leaf labels and written to a new tree file.
def addTaxonomyToTree(treeFileName, outFileName, IDfile = "/ncbiToTaxonomy.txt"):
    treeFile = open(treeFileName)
    tree = treeFile.readlines()[0].strip()
    treeFile.close()

    names = getLeafLabels(tree)

    NCBItoTax = {}

    IDlist = open(IDfile)
    IDlines = IDlist.readlines()
    for line in IDlines:
        l = line.strip()
        parts = l.split()
        #print(parts)
        if len(parts) > 1:
            NCBItoTax[parts[0]] = '_'.join(parts[1:])
            #print('b')
    IDlist.close()

    newIDs = {}
    taxonCounts = {}
    toAdd = {}

    ignore = False
    for name in names:
        if name.split('_')[0] in NCBItoTax:
            tax = NCBItoTax[name.split('_')[0]]
            newIDs[name] = name + '_' + tax
            if tax in taxonCounts:
                taxonCounts[tax] += 1
            else:
                taxonCounts[tax] = 1
        else:
            try:
                time.sleep(0.5)
                handle = Entrez.efetch(db = "protein", id = name.split('_')[0], rettype = "gb", retmode = "text")
                record = next(SeqIO.parse(handle, "genbank"))
                handle.close()
                tax = "_".join(record.annotations["taxonomy"] + [record.annotations["organism"].replace(' ', '_').replace('.', '')])
                if tax == "":
                    tax = "Bacteria_unclassified"
                newIDs[name] = name + '_' + tax
                if tax in taxonCounts:
                    taxonCounts[tax] += 1
                else:
                    taxonCounts[tax] = 1
                print("New lookup:",name, "->", tax)
                toAdd[name] = tax
            except:
                print("No record:", name)
                if not ignore:
                    ok = input("OK? (y/n/all): ")
                    if ok == 'y':
                        toAdd[name] = name
                    elif ok == 'all':
                        ignore = True
                        toAdd[name] = name
                if ignore:
                    toAdd[name] = name
                newIDs[name] = name
                if name in taxonCounts:
                    taxonCounts[name] += 1
                else:
                    taxonCounts[name] = 1

    if len(toAdd) > 0:
        IDlist = open("/ncbiToTaxonomy.txt", 'a')
        for ID in toAdd:
            line = ID + '\t' + toAdd[ID] + '\n'
            IDlist.write(line)
        IDlist.close()

    singles = []
    for taxon in taxonCounts:
        if taxonCounts[taxon] == 1:
            singles.append(taxon)
    #print(singles)
    for taxon in singles:
        taxonCounts.pop(taxon)
    for name in names:
        newName = newIDs[name]
        parts = newName.split('_')
        tax = '_'.join(parts[1:])
        if parts[0] in ["WP", "XP"]:
            tax = '_'.join(parts[2:])
        else:
            tax = '_'.join(parts[1:])
        #print(tax)
        if tax in taxonCounts:
            newName = newName + "_" + str(taxonCounts[tax])
            taxonCounts[tax] -= 1
            #print('a')
        tree = tree.replace(name, newName)

    out = open(outFileName, 'w')
    out.write(tree)
    out.close()

def writeLabelledTree(tree, outFileName, getTidFxn):
    linDict = {}
    nameDict = {}
    rankDict = {}
    treeLines = []
    taxLines = []
    taxTree = tree.copy()
    for node in taxTree.traverse():
        if node.is_leaf():
            parts = node.name.split('|')
            taxid = getTidFxn(node.name)
            gen_spe = nameDict.get(taxid)
            if not gen_spe:
                gen_spe = ncbi.get_taxid_translator([taxid])[int(taxid)]
            taxLine = '\t%s [&gen_spe="%s"' %(node.name, gen_spe)

            lin = linDict.get(taxid)
            if not lin:
                lin = ncbi.get_lineage(taxid)
                linDict[taxid] = lin
            notIn = [x for x in lin if x not in nameDict]
            names = ncbi.get_taxid_translator(notIn)
            ranks = ncbi.get_rank(notIn)
            for name in names:
                nameDict[name] = names[name]
                rankDict[name] = ranks[name]
            lineage = dict([(rankDict[x], nameDict[x]) for x in lin])
            for rank in ['class', 'phylum', 'order', 'family', 'kingdom', 'superkingdom']:
                if rank in lineage:
                    taxLine += ' tax_%s="%s"' %(rank, lineage[rank])
            taxLine += ' fullTax="%s"' %'_'.join([nameDict[x] for x in lin])
            taxLine += ']'
            taxLines.append(taxLine)
        else:
            if node.name:
                if '/' in node.name:
                    write_format = 1
                    aLRT, UFBoot = node.name.split('/')
                    node.name = '[&UFBoot=%.2f,aLRT=%.2f]' %(float(UFBoot), float(aLRT))
                else:
                    write_format = 0
                    node.support = float(node.name)

    newick_text = tree.write(format=write_format)
    if write_format:
        newick_text = re.sub('_&UFBoot_(\d+\.\d\d)_aLRT_(\d+\.\d\d)_', '[&UFBoot=\\1,aLRT=\\2]', newick_text)
    treeLines.append("\ttree %s = [&R] %s" %(outFileName, newick_text))

    out = open(outFileName, 'w')
    out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(taxLines))
    out.write('\n'.join(taxLines))
    out.write(';\nend;\n')
    out.write('begin trees;\n%s\nend;' %'\n'.join(treeLines))
    out.close()
