'''
analyse CDS with respect to structural mappings at protein level
'''
import argparse
import csv
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from joblib import Parallel, delayed

def makeCodonTable():
    '''
    make codon table
    '''
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    return codon_table

def mut(codon, pos, base):
    '''
    given a codon mutate position "pos" to base
    '''
    codon = list(codon)
    codon[pos] = base
    return ''.join(codon)

def nNnS(cds, codon_table):
    '''
    given a sequence calculate the number of possible nonsynomyous sites and synonymous sites
    '''
    bases = "TCAG"
    nN = 0
    nS = 0
    for AApos, codon in cds.iteritems():
        if codon in codon_table.keys():
            AA = codon_table[codon]
            for pos in range(len(codon)):
                codons = map(lambda base: codon_table[mut(codon, pos, base)], [b for b in bases if b != codon[pos]])
                nS += sum([ c == AA for c in codons])
                nN += sum([ c != AA for c in codons ])
    return {'nN': nN, 'nS': nS}
            
def splitCDSintoCodes(cds, splitlen = 3):
    '''
    split supplied CDS into group of 3 (i.e. code) for mapping to residues
    '''
    return {1 + (i+1)/3: cds[i:i+splitlen] for i in range(0, len(cds), splitlen) }
   
def findbounds(listOfData, category):
    '''
    find bounds of continuous segments of listOfData == category. Return pairs of indices of upper and lower bound of each segment.
    '''
    output = []
    cont = False
    start = None
    end = None
    for i in range(len(listOfData)):
        if listOfData[i] == category:
            if cont is False and start is None:
                start = i + 1
            cont = True
        else:    
            if cont is True:
                end = i
                if start != end:
                    output.append([start, end])
                else:
                    output.append(start)
                start = None
                cont = False
    if end is None and start is not None:
       output.append([start, len(listOfData)])
    return output    

def fetchFlank(cdsGrouped, AApos, flank = 1, before = True):
    '''
    fetch the CDS sequence before or after a AApos, length specified by "flank"
    '''
    if before is True and AApos > 1:
        sequence = cdsGrouped[AApos - 1]
        return sequence[-flank:]
        #return sequence[len(sequence) - 1]
    if before is True and AApos == 1:
        return ''
    if before is False and AApos < cdsGrouped.keys()[-1]:
        sequence = cdsGrouped[AApos + 1]
        return sequence[:flank]
        #return sequence[0]
    if before is False and AApos == cdsGrouped.keys()[-1]:
        return ''

def groupCDSbyRegion(cdsGrouped, regions):
    '''
    split grouped CDS (from function "splitCDSintoCodes" by regional mapping
    '''
    regionLevels = list(set(regions))
    regionBounds = {level: findbounds(regions, level) for level in regionLevels}
    output = {k: [] for k in regionLevels}
    try:
        for level in regionLevels:
            for bound in regionBounds[level]:
                if bound.__class__ is int:
                    minus1 = fetchFlank( cdsGrouped, bound, before = True )
                    plus1 = fetchFlank( cdsGrouped, bound, before = False )
                    output[level].append( "".join([minus1, cdsGrouped[bound], plus1]) )
                if bound.__class__ is list:
                    minus1 = fetchFlank( cdsGrouped, bound[0], before = True )
                    plus1 = fetchFlank( cdsGrouped, bound[1], before = False )
                    for j in range(bound[0], bound[1] + 1):
                        output[level].append( "".join([fetchFlank( cdsGrouped, j, before = True ), cdsGrouped[j], fetchFlank( cdsGrouped, j, before = False) ]) )
    except KeyError:
        print("KeyError")
        output = None

    return output

def countContexts(cdsSet, group, context):
    '''
    further process output of groupCDSbyRegion. count the number of occurences of context in cdsSet[category == group].
    '''
    if group in cdsSet.keys():
        sequences = cdsSet[group]
        output = [context in sequence for sequence in sequences]
    else:
        output = []
    return output    

def getEntry(seqset, entryID):
    '''
    return sequence from seq_entry
    '''
    return seqset[entryID].seq

def zoomvarCDSmap(uniprot, uniprotCDS, signatures, outfile, codon_table):
    '''
    overall pipeline, pass to parallelisation
    '''
    if '|' in uniprot:
        uniprotID = uniprot.split('|')[1]
    else:
        uniprotID = uniprot
    print(uniprotID)
    cds = splitCDSintoCodes(str( getEntry(uniprotCDS, uniprot) ))
    nNnSdict = nNnS(cds, codon_table)
    len_tot = len(cds)
    cdsAll = groupCDSbyRegion(cds, ["yes"] * len_tot)
    if type(cdsAll) != type(None):
        tot_sig_bool = {context[0]: [countContexts(cdsAll, "yes", c) for c in context] for context in signatures}
        tot_sig_bool = {context[0]: [any([ tot_sig_bool[context[0]][j][i] for j in range(len(tot_sig_bool[context[0]]))]) for i in range(len(tot_sig_bool[context[0]][0]))] for context in signatures}
        tot_sig = {context[0]: sum(tot_sig_bool[context[0]]) for context in signatures}
        out = [uniprotID, nNnSdict['nN'], nNnSdict['nS'], len_tot]
        out = out + [ tot_sig[context[0]] for context in signatures ]
        with open(outfile, 'a') as outstream:
            wr = csv.writer(outstream)
            wr.writerow(out)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Script for counting occurrences of each of the 32 SBS contexts in Uniprot Coding sequence set")
    parser.add_argument('cds', help = 'CDS fasta file')
    parser.add_argument('outfile', help = 'desired output filename (full path)')
    parser.add_argument('-ncore', '--num_core', type = int, default = 16,\
                        help = 'number of cores for parallelisation')
    args = parser.parse_args()

    codon_table = makeCodonTable()
    signatures = [['ACA','TGT'],['ACC','GGT'],['ACG','CGT'],['ACT','AGT'],['CCA','TGG'],['CCC','GGG'],['CCG','CGG'],['CCT','AGG'],['GCA','TGC'],['GCC','GGC'],['GCG','CGC'],['GCT','AGC'],['TCA','TGA'],['TCC','GGA'],['TCG','CGA'],['TCT','AGA'],['ATA','TAT'],['ATC','GAT'],['ATG','CAT'],['ATT','AAT'],['CTA','TAG'],['CTC','GAG'],['CTG','CAG'],['CTT','AAG'],['GTA','TAC'],['GTC','GAC'],['GTG','CAC'],['GTT','AAC'],['TTA','TAA'],['TTC','GAA'],['TTG','CAA'],['TTT','AAA']]
    uniprotCDS = SeqIO.index(args.cds, "fasta")
    uniprotIDs = [entry for entry in uniprotCDS]
    out = ['Uniprot', 'nN', 'nS', 'len_tot']
    with open(args.outfile, 'w') as outstream:
        wr = csv.writer(outstream)
        wr.writerow( out + [i[0] for i in signatures] )

    Parallel(n_jobs=args.num_core)(delayed(zoomvarCDSmap)(uniprot, uniprotCDS, \
        signatures, args.outfile, codon_table) for uniprot in uniprotIDs)        
