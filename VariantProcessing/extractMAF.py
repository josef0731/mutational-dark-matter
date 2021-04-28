'''
analyse on a per-sample basis the enrichment of APOBEC3-signature mutations at protein region level
'''
import argparse
import csv

def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]

def extractMafCols(maf, cols):
    '''
    straightforward extraction of the columns indicated from the MAF table. Assumed header (i.e. first row of maf) contains the column names.
    '''
    colsToExtract = [i for i in range(len(maf[0])) if maf[0][i] in cols]
    return [[entry[j] for j in colsToExtract] for entry in maf]

def tidyProteinsList(proteinsList, maf):
    '''
    tidy uniprot IDs from MAF file by trying to use other mappings from other columns
    '''
    for i in range(len(proteinsList)):
        if proteinsList[i] == '':
            uniprot = maf[ i + 1 ][[entry for entry in range(len(maf[0])) if maf[0][entry] == 'i_HGNC_UniProtIDsuppliedbyUniProt'][0]]
            if uniprot == '':
                uniprot = maf[ i + 1 ][[entry for entry in range(len(maf[0])) if maf[0][entry] == 'i_dbNSFP_Uniprot_acc'][0]]
                uniprot = uniprot.split(';')
                if len(uniprot) > 1:
                    uniprot_names = maf[ i + 1 ][[entry for entry in range(len(maf[0])) if maf[0][entry] == 'i_dbNSFP_Uniprot_acc'][0] + 1]
                    uniprot_names = uniprot_names.split(';')
                    for j in range(len(uniprot_names)-1, -1, -1):
                        if uniprot_names[j] == ".":
                            del uniprot[j]
                    if len(uniprot) == 0:
                        uniprot = ''
                else:
                    proteinsList[i] = uniprot
            else:                
                proteinsList[i] = uniprot       
    return proteinsList

def tidySample(maf, countsTb, \
    proteinColKey = 'SwissProt_acc_Id', positionColKey = 'UniProt_AApos'):
    '''
    tidy up the maf
    - tidy up uniprot IDs (call tidyProteinsList)
    - remove protein not mapped in zoomvar
    - remove protein &residues without regional mapping
    output:
    1) the cleaned maf
    2) list of unique Uniprots for the affected proteins
    '''
        
    proteinCol = [entry for entry in range(len(maf[0])) if maf[0][entry] == proteinColKey][0]
    positionCol = [entry for entry in range(len(maf[0])) if maf[0][entry] == positionColKey][0]
    proteinsList = [entry[proteinCol] for entry in maf[1:]]
    proteinsList = tidyProteinsList(proteinsList, maf)
    
    #remove those entries whose proteins/positions cannot map to zoomvar / structural regions
    for i in range(len(proteinsList)-1, -1, -1):
        if proteinsList[i] == '' or ';' in proteinsList[i] or '-' in proteinsList[i] or proteinsList[i] not in countsTb.keys():
            del proteinsList[i]
            del maf[i+1]
        else:
            maf[i+1][proteinCol] = proteinsList[i]

    proteinsList = unique(proteinsList)
    return [maf, proteinsList]        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Subsetting columns out of MAF files. Written for use with TCGA MAFs from Broad GDSC")
    parser.add_argument('maf', help = 'mutation annotation file')
    parser.add_argument('countsTb', help = 'per-protein counts of signature mutations at different regions (.csv)')
    parser.add_argument('outfile', help = 'destination to output')
    parser.add_argument('columns', nargs='+', help='list of columns to extract from maf. input as space-separated, unquoted texts')
    args = parser.parse_args()

    AAcountsTb = []
    with open(args.countsTb, 'r') as instream:
        read = csv.reader(instream)
        for row in read:
            AAcountsTb.append(row)

    countsTb = {}
    for row in AAcountsTb[1:]:
        countsTb[row[0]] = {AAcountsTb[0][j]: row[j] for j in range(len(row))}
        
    del AAcountsTb    

    maf = []
    with open(args.maf, 'r') as instream:
        read = csv.reader(instream, delimiter="\t")
        for row in read:
            maf.append(row)
    maf_cleaned = extractMafCols(maf, args.columns)
    maf_cleaned2 = tidySample(maf_cleaned, countsTb)[0]     
    with open(args.outfile, 'w') as outstream:
        wr = csv.writer(outstream, delimiter="\t")
        for row in maf_cleaned2:
            wr.writerow(row)
        
