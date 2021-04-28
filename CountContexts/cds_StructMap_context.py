'''
analyse CDS with respect to structural mappings at protein level
'''
import argparse
import csv
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
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
    return {1 + (i)/3: cds[i:i+splitlen] for i in range(0, len(cds), splitlen) }
   
def classifyRegion(regionMap, statCol, cutoff, below = 'core', above = 'surface'):
    '''
    classify region according to a statistic with supplied cutoff
    '''
    output = []
    for residue in regionMap:
        if residue[statCol] != 'None':
            if float(residue[statCol]) < cutoff:
                output.append(below)
            else:
                output.append(above)
        else:
            output.append('None')
    return output    
    
def classifyInteract(regions, regionMap, interactCol):
    '''
    classify whether residues are of type "interact". If yes, output an updated (to "interact") version of "regionCol"
    '''
    interact = classifyRegion(regionMap, statCol = interactCol, cutoff = 1, below = "No", above = "interact")
    for o in range(len(interact)):
        if interact[o] == "interact":
            regions[o] = "interact"
    return regions        

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
        print "KeyError"
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

#def zoomvarCDSmap(uniprot, uniprotCDS, zoomvar, disCol, qsasaCol, qsasa_cutoff, interaction, tcw, cytosine, outfile, codon_table):
def zoomvarCDSmap(uniprot, uniprotCDS, zoomvar, disCol, qsasaCol, qsasa_cutoff, interaction, tcw, outfile, codon_table):
    '''
    overall pipeline, pass to parallelisation
    '''
    uniprotID = uniprot.split("|")[1]
    print uniprotID
    cds = splitCDSintoCodes(str( getEntry(uniprotCDS, uniprot) ))
    nNnSdict = nNnS(cds, codon_table)
    zoomvarPart = [entry for entry in zoomvar if entry[1] == uniprotID]
    structure = classifyRegion(zoomvarPart, qsasaCol, qsasa_cutoff)
    structure = classifyInteract(structure, zoomvarPart, interaction)
    cdsByStructure = groupCDSbyRegion(cds, structure)
    len_tot = len(structure)
    cdsAll = groupCDSbyRegion(cds, ["yes"] * len_tot)
    if type(cdsAll) != type(None):
        tot_tcw_bool = [countContexts(cdsAll, "yes", context) for context in tcw]
        tot_tcw = sum([any([tot_tcw_bool[j][i] for j in range(len(tot_tcw_bool))]) for i in range(len(tot_tcw_bool[0]))])
        #tot_c_bool = [countContexts(cdsAll, "yes", context) for context in cytosine]
        #tot_c = sum([any([tot_c_bool[j][i] for j in range(len(tot_c_bool))]) for i in range(len(tot_c_bool[0]))])

    if type(cdsByStructure) != type(None):
        len_surf = structure.count("surface")
        len_core = structure.count("core")
        len_interact = structure.count("interact")
        core_tcw_bool = [countContexts(cdsByStructure, "core", context) for context in tcw]
        core_tcw = sum([any([core_tcw_bool[j][i] for j in range(len(core_tcw_bool))]) for i in range(len(core_tcw_bool[0]))])
        #core_c_bool = [countContexts(cdsByStructure, "core", context) for context in cytosine]
        #core_c = sum([any([core_c_bool[j][i] for j in range(len(core_c_bool))]) for i in range(len(core_c_bool[0]))])
        surf_tcw_bool = [countContexts(cdsByStructure, "surface", context) for context in tcw]
        surf_tcw = sum([any([surf_tcw_bool[j][i] for j in range(len(surf_tcw_bool))]) for i in range(len(surf_tcw_bool[0]))])
        #surf_c_bool = [countContexts(cdsByStructure, "surface", context) for context in cytosine]
        #surf_c = sum([any([surf_c_bool[j][i] for j in range(len(surf_c_bool))]) for i in range(len(surf_c_bool[0]))])
        interact_tcw_bool = [countContexts(cdsByStructure, "interact", context) for context in tcw]
        interact_tcw = sum([any([interact_tcw_bool[j][i] for j in range(len(interact_tcw_bool))]) for i in range(len(interact_tcw_bool[0]))])
        #interact_c_bool = [countContexts(cdsByStructure, "interact", context) for context in cytosine]
        #interact_c = sum([any([interact_c_bool[j][i] for j in range(len(interact_c_bool))]) for i in range(len(interact_c_bool[0]))])

    disorder = [residue[disCol] for residue in zoomvarPart]
    cdsByDisorder = groupCDSbyRegion(cds, disorder)
    if type(cdsByDisorder) != type(None):
        len_intradis = disorder.count('intra_dis')
        len_interdis = disorder.count('inter_dis')
        len_intraord = disorder.count('intra_ord')
        len_interord = disorder.count('inter_ord')
        intradis_tcw_bool = [countContexts(cdsByDisorder, "intra_dis", context) for context in tcw]
        intradis_tcw = sum([any([intradis_tcw_bool[j][i] for j in range(len(intradis_tcw_bool))]) for i in range(len(intradis_tcw_bool[0]))])
        #intradis_c_bool = [countContexts(cdsByDisorder, "intra_dis", context) for context in cytosine]
        #intradis_c = sum([any([intradis_c_bool[j][i] for j in range(len(intradis_c_bool))]) for i in range(len(intradis_c_bool[0]))])
        interdis_tcw_bool = [countContexts(cdsByDisorder, "inter_dis", context) for context in tcw]
        interdis_tcw = sum([any([interdis_tcw_bool[j][i] for j in range(len(interdis_tcw_bool))]) for i in range(len(interdis_tcw_bool[0]))])
        #interdis_c_bool = [countContexts(cdsByDisorder, "inter_dis", context) for context in cytosine]
        #interdis_c = sum([any([interdis_c_bool[j][i] for j in range(len(interdis_c_bool))]) for i in range(len(interdis_c_bool[0]))])
        intraord_tcw_bool = [countContexts(cdsByDisorder, "intra_ord", context) for context in tcw]
        intraord_tcw = sum([any([intraord_tcw_bool[j][i] for j in range(len(intraord_tcw_bool))]) for i in range(len(intraord_tcw_bool[0]))])
        #intraord_c_bool = [countContexts(cdsByDisorder, "intra_ord", context) for context in cytosine]
        #intraord_c = sum([any([intraord_c_bool[j][i] for j in range(len(intraord_c_bool))]) for i in range(len(intraord_c_bool[0]))])
        interord_tcw_bool = [countContexts(cdsByDisorder, "inter_ord", context) for context in tcw]
        interord_tcw = sum([any([interord_tcw_bool[j][i] for j in range(len(interord_tcw_bool))]) for i in range(len(interord_tcw_bool[0]))])
        #interord_c_bool = [countContexts(cdsByDisorder, "inter_ord", context) for context in cytosine]
        #interord_c = sum([any([interord_c_bool[j][i] for j in range(len(interord_c_bool))]) for i in range(len(interord_c_bool[0]))])

    if type(cdsAll) == type(None) or type(cdsByDisorder) == type(None) or type(cdsByStructure) == type(None):
        with open('failedCDS_proteinLevel_NEW2.csv', 'a') as outstream:
            wr = csv.writer(outstream)
            wr.writerow([uniprotID])
    else:
        with open(outfile, 'a') as outstream:
            wr = csv.writer(outstream)
            wr.writerow([uniprotID, nNnSdict['nN'], nNnSdict['nS'], len_tot, len_surf, len_core, len_interact, \
                len_intradis, len_interdis, len_intraord, len_interord, tot_tcw, \
		#tot_c
                core_tcw, surf_tcw, interact_tcw, intradis_tcw, \
		#intradis_c,core_tcw, core_c, surf_tcw, surf_c, interact_tcw, interact_c, intradis_tcw, \#intradis_c
                interdis_tcw, intraord_tcw, interord_tcw])
#                interdis_tcw, interdis_c, intraord_tcw, intraord_c, interord_tcw, interord_c])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Script for the annotation of Uniprot Coding sequence by Zoomvar structural mapping and APOBEC sequence context counts")
    parser.add_argument('zoomvar', help = 'Zoomvar mapping tsv')
    parser.add_argument('cds', help = 'CDS fasta file')
    parser.add_argument('outfile', help = 'desired output filename (full path)')
    parser.add_argument('-sig', '--sig', nargs="+", help = 'context to count')
    parser.add_argument('-ncore', '--num_core', type = int, default = 16,\
                        help = 'number of cores for parallelisation')
    parser.add_argument('-qsasa', '--qsasa', type = float, default = 0.15,\
                        help = 'qsasa cut-off to define surface/core')
    args = parser.parse_args()

    codon_table = makeCodonTable()

    zoomvar = []
    with open(args.zoomvar, 'r') as instream:
        read = csv.reader(instream, delimiter="\t")
        for row in read:
            zoomvar.append(row)

    for i in range(len(zoomvar[0])):
        if zoomvar[0][i] == "Q(SASA)":
            qsasaCol = i
        if zoomvar[0][i] == "position":
            residue = i
        if zoomvar[0][i] == "disorder_cat":
            disCol = i
        if zoomvar[0][i] == "num_interactions":
            interaction = i
    
    #tcw = ['TCA', 'TCT', 'TGA', 'AGA']
    tcw = args.sig
    #cytosine = ['C', 'G']
    uniprotCDS = SeqIO.index(args.cds, "fasta")
    uniprotIDs = [entry for entry in uniprotCDS]
    with open(args.outfile, 'w') as outstream:
        wr = csv.writer(outstream)
        wr.writerow(['Uniprot', 'nN', 'nS', 'len_tot', 'len_surface', 'len_core', 'len_interact', 'len_intradis', 'len_interdis', 'len_intraord', 'len_interord', 'tot_sig', 'core_sig', 'surface_sig', 'interact_sig', 'interdis_sig', 'intradis_sig', 'interord_sig', 'intraord_sig'])

#        wr.writerow(['Uniprot', 'nN', 'nS', 'len_tot', 'len_surface', 'len_core', 'len_interact', 'len_intradis', 'len_interdis', 'len_intraord', 'len_interord', 'tot_tcw', 'tot_c', 'core_tcw', 'core_c', 'surface_tcw', 'surface_c', 'interact_tcw', 'interact_c', 'interdis_tcw', 'interdis_c', 'intradis_tcw', 'intradis_c', 'interord_tcw', 'interord_c', 'intraord_tcw', 'intraord_c'])
    #Parallel(n_jobs=args.num_core)(delayed(zoomvarCDSmap)(uniprot, uniprotCDS, zoomvar, disCol, qsasaCol, args.qsasa, interaction, tcw, cytosine, args.outfile, codon_table) for uniprot in uniprotIDs)        
    Parallel(n_jobs=args.num_core)(delayed(zoomvarCDSmap)(uniprot, uniprotCDS, zoomvar, disCol, qsasaCol, args.qsasa, \
        interaction, tcw, args.outfile, codon_table) for uniprot in uniprotIDs)        
