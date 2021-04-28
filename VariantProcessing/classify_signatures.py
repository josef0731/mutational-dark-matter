# classify mutations into
# - UV
# - 5-FU
# - Platinum
import argparse
import csv
import re

def isApobecSig(context, ref, t1, t2):
    """
    Tell whether the mutation is in APOBEC signature
    context - the 20-bp 5' and 3' of the mutation
              (i.e. 41bp sequence, mutated position at the 21st position.)
    ref - reference (normal) allele
    t1  - tumour allele 1
    t2  - tumour allele 2
    """
    mut_pos = 20 # 0-based
    if (context[(mut_pos-1):(mut_pos+2)].upper() in ["TCT", "TCA", "AGA", "TGA"]) and \
        ((ref == "C" and (t1 == "T" or t2 == "T" or t1 == "G" or t2 == "G")) or \
        (ref == "G" and (t1 == "A" or t2 == "A" or t1 == "C" or t2 == "C"))):
        return 1
    else:
        return 0

def isAgingSig(context, ref, t1, t2):
    """
    Tell whether the mutation is in Aging signature
    context - the 20-bp 5' and 3' of the mutation
              (i.e. 41bp sequence, mutated position at the 21st position.)
    ref - reference (normal) allele
    t1  - tumour allele 1
    t2  - tumour allele 2
    """
    mut_pos = 20 # 0-based
    if (context[mut_pos:(mut_pos+2)].upper() in ["CG", "GC"]) and \
        ((ref == "C" and (t1 == "T" or t2 == "T")) or (ref == "G" and (t1 == "A" or t2 == "A"))):
        return 1
    else:
        return 0

def is5FUSig(context, ref, t1, t2):
    """
    Tell whether the mutation is in 5-FU signature
    context - the 20-bp 5' and 3' of the mutation
              (i.e. 41bp sequence, mutated position at the 21st position.)
    ref - reference (normal) allele
    t1  - tumour allele 1
    t2  - tumour allele 2
    """
    mut_pos = 20 # 0-based
    if (context[(mut_pos-1):(mut_pos+2)].upper() in ["CTT", "AAG"]) and \
        ((ref == "T" and (t1 == "G" or t2 == "G")) or (ref == "A" and (t1 == "C" or t2 == "C"))):
        return 1
    else:
        return 0

def isPOLESig(context, ref, t1, t2):
    """
    Tell whether the mutation is in POLE signature
    context - the 20-bp 5' and 3' of the mutation
              (i.e. 41bp sequence, mutated position at the 21st position.)
    ref - reference (normal) allele
    t1  - tumour allele 1
    t2  - tumour allele 2
    """
    mut_pos = 20 # 0-based
    if (context[(mut_pos-1):(mut_pos+2)].upper() in ["TCT", "AGA"]) and \
        ((ref == "C" and (t1 == "A" or t2 == "A")) or (ref == "G" and (t1 == "T" or t2 == "T"))):
        return 1
    elif (context[(mut_pos-1):(mut_pos+2)].upper() in ["TCG", "CGA"]) and \
        ((ref == "C" and (t1 == "T" or t2 == "T")) or (ref == "G" and (t1 == "A" or t2 == "A"))):
        return 1
    else:
        return 0

def isUVSig(ref, t1, t2):
    """
    Tell whether the mutation is in UV signature
    ref - reference (normal) allele. string of length 2 (since it is a doublet substitution)
    t1  - tumour allele 1. string of length 2 (since it is a doublet substitution)
    t2  - tumour allele 2. string of length 2 (since it is a doublet substitution)
    """
    if ref == "CC" and (t1 == "TT" or t2 == "TT"):
        return 1
    elif ref == "GG" and (t1 == "AA" or t2 == "AA"):
        return 1
    else:
        return 0

def isPlatSig(ref, t1, t2):
    """
    Tell whether the mutation is in Platinum signature
    ref - reference (normal) allele. string of length 2 (since it is a doublet substitution)
    t1  - tumour allele 1. string of length 2 (since it is a doublet substitution)
    t2  - tumour allele 2. string of length 2 (since it is a doublet substitution)
    """
    if ref == "CT" and (t1 == "AA" or t2 == "AA"):
        return 1
    elif ref == "CT" and (t1 == "AC" or t2 == "AC"):
        return 1
    elif ref == "AG" and (t1 == "TT" or t2 == "TT"):
        return 1
    elif ref == "AG" and (t1 == "GT" or t2 == "GT"):
        return 1
    elif ref == "TG" and (t1 == "GT" or t2 == "GT"):
        return 1
    elif ref == "TC" and (t1 == "AA" or t2 == "AA"):
        return 1
    elif ref == "CA" and (t1 == "AC" or t2 == "AC"):
        return 1
    elif ref == "GA" and (t1 == "TT" or t2 == "TT"):
        return 1
    else:
        return 0

def isPlatSigSBS(context, ref, t1, t2):
    """
    Tell whether the mutation is in Platinum SBS signature (COSMIC signature 31)
    context - the 20-bp 5' and 3' of the mutation
              (i.e. 41bp sequence, mutated position at the 21st position.)
    ref - reference (normal) allele
    t1  - tumour allele 1
    t2  - tumour allele 2
    """
    mut_pos = 20 # 0-based
    if (context[(mut_pos-1):(mut_pos+2)].upper() in ["CCC", "GGG"]) and \
        ((ref == "C" and (t1 == "T" or t2 == "T")) or (ref == "G" and (t1 == "A" or t2 == "A"))):
        return 1
    elif (context[(mut_pos-1):(mut_pos+2)].upper() in ["CCT", "AGG"]) and \
        ((ref == "C" and (t1 == "T" or t2 == "T")) or (ref == "G" and (t1 == "A" or t2 == "A"))):
        return 1
    else:
        return 0

def getDoublet(refA, refB, t1A, t1B, t2A, t2B):
    """
    parse two separate SNVs into a doublet
    """
    return {'ref': refA + refB, 't1': t1A + t1B, 't2': t2A + t2B}

def makeCodonTable():
    '''
    make codon table
    '''
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))
    return codon_table

def parseDoubletConsqOnProt(protein_changes, structs, disords, origs, muts, varclasses, \
    context, codonTable, ref = None, t1 = None, t2 = None):
    """
    parse the consequence of doublet mutation on protein
    """
    protein_pos = [re.sub('[A-Za-z]', '', j.replace('.', '').replace('*', '')) for j in protein_changes]
    protein_pos = list(set(protein_pos))
    # first rewrite context so that both bases are highlighted as capital letters
    context_list = []
    for i in range(len(context)):
        if i == 21:
            context_list.append(context[i].upper())
        else:
            context_list.append(context[i])
    context = ''.join(context_list)
    # test if affecting 1 or 2 protein positions
    if len(protein_pos) == 2:
        # no need to re-evaluate mutated AA. But need to parse
        # mapping to struct and disord regions
        # if one residue map to a region but the other is None,
        # remove the None. If mapped to two regions, return both regions.
        structs = list(set(structs))
        if len(structs) == 1 and structs[0] == 'None':
            structs = ['None']
        elif len(structs) == 2 and 'None' in structs:
            structs = [i for i in structs if i is not 'None']
        disords = list(set(disords))
        if len(disords) == 1 and disords[0] == 'None':
            disords = ['None']
        elif len(disords) == 2 and 'None' in disords:
            disords = [i for i in disords if i is not 'None']
        if varclasses[0] == varclasses[1]:
            varclass = varclasses[0]
        else:
            if 'Nonsense_Mutation' in varclasses:
                varclass = 'Nonsense_Mutation'
            elif sum([i == 'Silent' for i in varclasses]) == 1:
                varclasses.remove('Silent')
                varclass = ';'.join(varclasses)
            else:
                varclass = ';'.join(varclasses)
        if origs[0] == origs[1]:
            orig = origs[0]
        else:
            orig = ''.join(origs)
        if muts[0] == muts[1]:
            mut = muts[0]
        else:
            mut = ''.join(muts)
        return { 'Protein_Change': ';'.join(protein_changes), 'context': context, \
            'structs': ';'.join(structs), 'disords': ';'.join(disords), 'varclass': varclass, \
            'origs': orig, 'muts': mut}
    elif len(protein_pos) == 1:
        # affecting same AA position but now has 2 different mut AA. Need
        # to parse actual mutated AA but no need to parse struct and disord mapping
        is_substitution = sum([ protein_change.startswith('p.') for protein_change in protein_changes ])
        if is_substitution == 0:
            # not substitution. in this case Protein_Change is None. And struct and disord is None.
            if len(list(set(origs))) == 1 and origs[0] == '':
                origs = ['.']
            if len(list(set(muts))) == 1 and muts[0] == '':
                muts = ['.']
            return { 'Protein_Change': '.', 'context': context, 'origs': ''.join(origs), 'muts': ''.join(muts), \
                'structs': 'None', 'disords': 'None', 'varclass': ';'.join(list(set(varclasses))) }
        protein_changes = protein_changes[0]
        origAA = protein_changes[2]
        wt = context[19:22].upper()
        if codonTable[wt] == origAA:
            if ref == t1:
                mut = context[19] + t2
            else:
                mut = context[19] + t1
        else:
            wt = context[20:23].upper()
            if ref == t1:
                mut = context[20] + t2
            else:
                mut = context[20] + t1
        mut = mut.upper()
        mutAA = codonTable[mut]
        if mutAA == "*":
            varclass = 'Nonsense_Mutation'
        elif origAA == mutAA:
            varclass = 'Silent'
        elif origAA != mutAA:
            varclass = 'Missense_Mutation'
        else:
            varclass = ';'.join(varclasses)
        return { 'Protein_Change': 'p.' + origAA + protein_pos[0] + mutAA, 'context': context, \
            'structs': structs[0], 'disords': disords[0], 'varclass': varclass, \
            'origs': origAA, 'muts': mutAA }
#    else:
#        return { 'Protein_Change': ';'.join(protein_changes), 'context': context, \
#            'structs': ';'.join(structs), 'disords': ';'.join(disords), 'varclass': varclass, \
#            'origs': orig, 'muts': mut }

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This script read MAF files and annotate variants for selected mutational signatures (APOBEC, UV, POLE, Aging, 5-FU, Platinum)")
    parser.add_argument('--mafIn', help = 'maf input (TCGA firehose MAF standard format with annotation of sequence context around each mutated position)', type = str, required = True)
    parser.add_argument('--mafOut', help = 'maf output with annotation added', type = str, required = True)
    args = parser.parse_args()

    maf = []
    header = []
    with open(args.mafIn, 'r') as instream:
        read = csv.reader(instream, delimiter="\t")
        for row in read:
            if row[0] == 'Hugo_Symbol':
                for column in range(len(row)):
                    if row[column] == 'Chromosome':
                        chrCol = column
                    elif row[column] == "Start_Position":
                        startCol = column
                    elif row[column].startswith("End_"):
                        endCol = column
                    elif row[column] == "Strand":
                        strandCol = column
                    elif row[column] == "Reference_Allele":
                        refCol = column
                    elif row[column] == "Tumor_Seq_Allele1":
                        t1Col = column
                    elif row[column] == "Tumor_Seq_Allele2":
                        t2Col = column
                    elif row[column] == "Tumor_Sample_Barcode":
                        sampleCol = column
                    elif row[column] == "Protein_Change":
                        proteinCol = column
                    elif row[column] == "struct_cat":
                        structCol = column
                    elif row[column] == "disorder_cat":
                        disordCol = column
                    elif row[column] == "Variant_Classification":
                        varclassCol = column
                    elif row[column] == "origAA":
                        origCol = column
                    elif row[column] == "mutAA":
                        mutCol = column
                    elif row[column].startswith("CONTEXT"):
                        contextCol = column
                header = row
            else:
                maf.append(row)

    print "Read maf file " + args.mafIn
    codon_table = makeCodonTable()
    with open(args.mafOut, 'w') as outstream:
        wr = csv.writer(outstream, delimiter="\t")
        wr.writerow(header + ['ApobecSig', 'AgingSig', 'POLESig', 'UVSig', 'FiveFUSig', 'PlatDBS', 'PlatSBS', 'PlatSig'])
        chr = start = end = strand = ref = t1 = t2 = sample = None
        for row in maf:
            is_doublet = False
            is_consecutive = False
            if row[refCol] in ['A', 'T', 'C', 'G'] and row[t1Col] in ['A', 'T', 'C', 'G'] and \
                row[t2Col] in ['A', 'T', 'C', 'G']:
                FiveFUSig = is5FUSig(row[contextCol], row[refCol], row[t1Col], row[t2Col])
                POLESig = isPOLESig(row[contextCol], row[refCol], row[t1Col], row[t2Col])
                AgingSig = isAgingSig(row[contextCol], row[refCol], row[t1Col], row[t2Col])
                ApobecSig = isApobecSig(row[contextCol], row[refCol], row[t1Col], row[t2Col])
                PlatSigSBS = isPlatSigSBS(row[contextCol], row[refCol], row[t1Col], row[t2Col])
                if PlatSigSBS == 1:
                    PlatSig = 1
                else:
                    PlatSig = 0
                PlatSigDBS, UVSig = 'NA', 'NA'
            else:
                FiveFUSig, POLESig, AgingSig, ApobecSig, PlatSigSBS = 'NA', 'NA', 'NA', 'NA', 'NA'
                if len(re.sub('-', '', row[refCol])) == 2 and \
                    ( len(re.sub('-', '', row[t1Col])) == 2 and len(re.sub('-', '', row[t2Col])) == 2 ):
                    PlatSigDBS = isPlatSig(row[refCol], row[t1Col], row[t2Col])
                    if PlatSigDBS == 1:
                        PlatSig = 1
                    else:
                        PlatSig = 0
                    UVSig = isUVSig(row[refCol], row[t1Col], row[t2Col])
                    is_doublet = True
                else:
                    PlatSigDBS, UVSig, PlatSig = 'NA', 'NA', 'NA'
            if any([i is None for i in [chr, start, end, strand, ref, t1, t2, sample]]):
                PlatSigDBS, UVSig, PlatSig = 'NA', 'NA', 'NA'
            elif row[sampleCol] == sample and row[chrCol] == chr and row[strandCol] == strand and \
                start == end and (int(row[startCol]) == end + 1):
                if len(re.sub('-', '', row[refCol])) == 1 and len(re.sub('-', '', ref)) == 1 and \
                    len(re.sub('-', '', row[t1Col])) == 1 and len(re.sub('-', '', t1)) == 1 and \
                    len(re.sub('-', '', row[t2Col])) == 1 and len(re.sub('-', '', t2)) == 1:
                    is_consecutive = True
                    if previous_row is None and just_in_case is not None:
                        # this is actually a string of substitution longer than doublet
                        # in this case just output as separate base substitutions
                        wr.writerow(just_in_case)
                        is_consecutive = False
                    else:
                        previous_row[endCol] = row[endCol]
                        DB = getDoublet(ref, row[refCol], t1, row[t1Col], t2, row[t2Col])
                        previous_row[refCol] = DB['ref']
                        previous_row[t1Col] = DB['t1']
                        previous_row[t2Col] = DB['t2']
                        PlatSigDBS = isPlatSig(DB['ref'], DB['t1'], DB['t2'])
                        if PlatSigDBS == 1:
                            PlatSig = 1
                        else:
                            PlatSig = 0
                        if isPlatSig(DB['ref'], DB['t1'], DB['t2']) == 1:
                            previous_row[len(previous_row)-3] = '1'
                        elif isPlatSig(DB['ref'], DB['t1'], DB['t2']) == 0:
                            previous_row[len(previous_row)-3] = '0'
                        UVSig = isUVSig(DB['ref'], DB['t1'], DB['t2'])
                        if isUVSig(DB['ref'], DB['t1'], DB['t2']) == 1:
                            previous_row[len(previous_row)-5] = '1'
                        elif isUVSig(DB['ref'], DB['t1'], DB['t2']) == 0:
                            previous_row[len(previous_row)-5] = '0'
                    # erase SBS as this is a doublet variant
                        previous_row[len(previous_row)-2] = 'NA'
                        previous_row[len(previous_row)-4] = 'NA'
                        previous_row[len(previous_row)-6] = 'NA'
                        previous_row[len(previous_row)-7] = 'NA'
                        previous_row[len(previous_row)-8] = 'NA'
                        # parse protein consequence
                        protConsq = parseDoubletConsqOnProt(protein_changes = [protein_change, row[proteinCol]], \
                            structs = [struct, row[structCol]], disords = [disord, row[disordCol]], \
                            context = context, codonTable = codon_table, origs = [orig, row[origCol]], \
                            muts = [mut, row[mutCol]], varclasses = [varclass, row[varclassCol]], \
                            ref = DB['ref'], t1 = DB['t1'], t2 = DB['t2'])
                        previous_row[proteinCol] = protConsq['Protein_Change']
                        previous_row[contextCol] = protConsq['context']
                        previous_row[structCol] = protConsq['structs']
                        previous_row[disordCol] = protConsq['disords']
                        previous_row[varclassCol] = protConsq['varclass']
                        previous_row[origCol] = protConsq['origs']
                        previous_row[mutCol] = protConsq['muts']
                        wr.writerow(previous_row)
            elif is_doublet:
                if previous_row is not None:
                    wr.writerow(previous_row)
            else:
                if previous_row is not None:
                    wr.writerow(previous_row)
            chr = row[chrCol]
            start = int(row[startCol])
            end = int(row[endCol])
            strand = row[strandCol]
            ref = row[refCol]
            t1 = row[t1Col]
            t2 = row[t2Col]
            sample = row[sampleCol]
            protein_change = row[proteinCol]
            struct = row[structCol]
            disord = row[disordCol]
            context = row[contextCol]
            varclass = row[varclassCol]
            orig = row[origCol]
            mut = row[mutCol]
            for item in range(len(row)):
                if row[item] == '':
                    row[item] = '.'
            if is_consecutive:
                previous_row = None
                just_in_case = row + [str(ApobecSig), str(AgingSig), str(POLESig), str(UVSig), str(FiveFUSig), str(PlatSigDBS), str(PlatSigSBS), str(PlatSig)]
            else:
                previous_row = row + [str(ApobecSig), str(AgingSig), str(POLESig), str(UVSig), str(FiveFUSig), str(PlatSigDBS), str(PlatSigSBS), str(PlatSig)]
                just_in_case = None
