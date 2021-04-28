'''
analyse on a per-sample basis the enrichment of signature mutations at protein region level
'''
import argparse
import collections
import csv
import cPickle
from joblib import Parallel, delayed
import numpy as np
import random
import re
import scipy.stats 
from scipy.stats import binom
amino_acids = ['C','S','T','P','A','G','N','D','E','Q','H','R','K','M','I','L','V','F','Y','R']

def binomCDF(listOfProteins, countsTb, category, regions, context, nmut_all, nmut_region): 
    '''
    calculate binomial cdf over a list of proteins
    :param context: either 'tcw' or 'c'
    :param nmut_all: the number of mutations in the context in all regions
    :param nmut_region: the number of mutations in the context in a given region to calculate
    '''
    if '_' in category:
        category = re.sub('_', '', category)
        regions = [re.sub('_', '', region) for region in regions]
    if nmut_region > 0:
        n = nmut_all
        x = nmut_region
        size_region = float(sum([int(countsTb[protein][category + '_' + context]) \
                                 for protein in listOfProteins])) 
        size_all = float(sum([sum([int(countsTb[protein][region + '_' + context]) \
                                   for region in regions]) \
                              for protein in listOfProteins]))
        p = size_region / size_all
        return binom.cdf(x, n, p)
    else:
        return None

def apobecStructuralEnrichment(listOfProteins, countsTb, region, regions, nmut_sign, nmut_all, \
               context = 'tcw'):
    '''
    calculate density of APOBEC-signature mutations in a given region relative all regions combined
    this is normalised by the # of residues mutatable in an APOBEC context (i.e. "tcw") in the given region
    and in all regions combined
    '''
    if '_' in region:
        region = re.sub('_', '', region)
        regions = [re.sub('_', '', r) for r in regions]
    if nmut_all > 0:
        size_sign_region = sum([int(countsTb[protein][region + '_' + context]) \
                        for protein in listOfProteins])
        size_sign_all = sum([sum([int(countsTb[protein][region + '_' + context]) for region in regions]) \
                        for protein in listOfProteins])
        if size_sign_region > 0 and size_sign_all > 0:
            return (nmut_sign / float(size_sign_region)) / (nmut_all / float(size_sign_all))
        else:
            return None
    else:
        return None


def simNull_struct(niter, pos_mut, countsTb):
    '''
    simulation of null: for each iteration number of mutations in each protein involved is kept constant, but the position is randomised
    
    pos_mut is a counter object of number of mutations localise to each protein (uniprot ID as keys).

    Return a list of len=niter, dictionaries of how many of the sampled positions are in all combinations of structural regions(surface/core/interface) & seq context (sign/base/no)
    
    '''        
    lengths = {}
    for protein in pos_mut.keys():
        lengths[protein] = {'total': int(countsTb[protein]['core_sig']) \
                                     + int(countsTb[protein]['surface_sig']) \
                                     + int(countsTb[protein]['interact_sig']), \
                            'core_sig': int(countsTb[protein]['core_sig']), \
                            'surface_sig': int(countsTb[protein]['surface_sig']), \
                            'interact_sig': int(countsTb[protein]['interact_sig'])}
    #for each iteration: for each protein draw x number(s) from the list 1:len_tot, where x is the len of pos_mut specific to the proteincoun
    out = []
    out = Parallel(n_jobs=4)(delayed(samplePosition_struct)(i, lengths, pos_mut) for i in range(niter))
    #for i in range(niter):
    #    if i in range(0, 10000, 100):
    #        print "iteration " + str(i) + '...'
    #    out.append(samplePosition_struct(lengths, pos_mut))
    return out

def samplePosition_struct(i, lengths, pos_mut, replacement = False):
    '''
    sample positions from protein. Default without replacement
    Return a list of len=niter, dictionaries of how many of the sampled positions are in all combinations of structural regions(surface/core/interface) & seq context (sign/base/no)
    '''
    if i in range(0, 10000, 5):    
        print "iteration " + str(i) + '...'
    o = {key: 0 for key in ['proteins', 'core_sig', 'surface_sig', 'interact_sig']}
    o['proteins'] = []
    s = {}
    for protein in lengths.keys():
        if pos_mut[protein] <= lengths[protein]['total']:
            sample_n = pos_mut[protein]
        else:
            sample_n = lengths[protein]['total']
        if lengths[protein]['total'] > 0:
            if replacement is True:
                s[protein] = random.choice(range(1, lengths[protein]['total'] + 1), sample_n)
            else:
                s[protein] = random.sample(range(1, lengths[protein]['total'] + 1), sample_n)
        
    for protein in s.keys():
        for i in s[protein]:
            o['proteins'].append(protein)
            core_tcw = lengths[protein]['core_sig']
            surface_tcw = lengths[protein]['surface_sig']
            interact_tcw = lengths[protein]['interact_sig']
            if i < core_tcw:
                o['core_sig'] += 1
            elif i < core_tcw + surface_tcw:
                o['surface_sig'] += 1
            elif i < core_tcw + surface_tcw + interact_tcw:
                o['interact_sig'] += 1
    return o            


def simNull_disorder(niter, pos_mut, countsTb):
    '''
    simulation of null: for each iteration number of mutations in each protein involved is kept constant, but the position is randomised
    
    pos_mut is a dict, each entry is the number of mutations in each protein in protein_mut, uniprot IDs as keys.

    Return a list of len=niter, dictionaries of how many of the sampled positions are in all combinations of structural regions(surface/core/interface) & seq context (sign/base/no)
    
    '''        
    lengths = {}
    for protein in pos_mut.keys():
        lengths[protein] = {'total': int(countsTb[protein]['intradis_tcw']) \
                                     + int(countsTb[protein]['interdis_tcw']) \
                                     + int(countsTb[protein]['intraord_tcw']) \
                                     + int(countsTb[protein]['interord_tcw']), \
                            'intradis_tcw': int(countsTb[protein]['intradis_tcw']), \
                            'interdis_tcw': int(countsTb[protein]['interdis_tcw']), \
                            'intraord_tcw': int(countsTb[protein]['intraord_tcw']), \
                            'interord_tcw': int(countsTb[protein]['interord_tcw']) } \
    #for each iteration: for each protein draw x number(s) from the list 1:len_tot, 
    #where x is the len of pos_mut specific to the proteincoun
    out = []
    for i in range(niter):
        out.append(samplePosition_disorder(lengths, pos_mut))
    return out    

def samplePosition_disorder(lengths, pos_mut, replacement = False):
    '''
    sample positions from protein. Default without replacement
    Return a list of len=niter, dictionaries of how many of the sampled positions are in all combinations of structural regions(surface/core/interface) & seq context (sign/base/no)
    '''
    o = {key: 0 for key in ['proteins', 'intradis_tcw', 'interdis_tcw', \
                            'intraord_tcw', 'interord_tcw']}
    o['proteins'] = []
    s = {}
    for protein in lengths.keys():
        if pos_mut[protein] <= lengths[protein]['total']:
            sample_n = pos_mut[protein]
        else:
            sample_n = lengths[protein]['total']
        if lengths[protein]['total'] > 0:
            if replacement is True:
                s[protein] = random.choice(range(1, lengths[protein]['total'] + 1), sample_n)
            else:
                s[protein] = random.sample(range(1, lengths[protein]['total']+ 1), sample_n)

    for protein in s.keys():
        for i in s[protein]:
            o['proteins'].append(protein)
            intradis_tcw = lengths[protein]['intradis_tcw']
            interdis_tcw = lengths[protein]['interdis_tcw']
            intraord_tcw = lengths[protein]['intraord_tcw']
            interord_tcw = lengths[protein]['interord_tcw']
               
            if i < intradis_tcw:
                o['intradis_tcw'] += 1
            elif i < intradis_tcw + interdis_tcw:
                o['interdis_tcw'] += 1
            elif i < intradis_tcw + interdis_tcw + intraord_tcw:
                o['intraord_tcw'] += 1
            elif i < intradis_tcw + interdis_tcw + intraord_tcw + interord_tcw:
                o['interord_tcw'] += 1
    return o            

def sample_with_replacement(i, maf, nsample):
    if i in range(0, 10000, 5):
        print 'iteration ' + str(i) + ' ...'
    return [random.choice(maf[1:]) for _ in range(nsample)]

def bootstrap(niter, nsample, maf):
    '''
    bootstrap population: for each iteration in the list of mutations of the sample resample with replacement the same number of mutations as those in the subcategory to be computed

    straightforward sampling from maf

    nsample is the number of mutations to sample in each iteration

    Return a list of len=niter, dictionaries of the proteins, the positions and the CDS sequences (both the three-letter codes and five-letters (codon + +1 and -1 base)

    '''
    out=Parallel(n_jobs=4)(delayed(sample_with_replacement)(i, maf, nsample) for i in range(niter))
    #out = []
    #for i in range(niter):
    #    random.seed(i)
    #    if i in range(0, 10000, 100):
    #        print "iteration " + str(i) + '...'
    #    s = [random.choice(maf[1:]) for _ in range(nsample)] #sample with replacement
    #    out.append(s)
    return out

def getSignatureCountsFromMaf(population, region, regions, mapCol, signCol, sigVal, baseCol, altCol, proteinCol):
    '''
    get signature counts from lists of MAFs
    return a list of lists indicating the number of signature mutations and number of -ve ctrl mutations, and the unique list of Uniprot IDs mutated
    '''
    out = []
    for item in population:
        in_region = [item[i] for i in range(len(item)) if item[i][mapCol] == region]
        s_sign = sum([ entry[signCol] in sigVal for entry in in_region] )
        s_sign_all = sum([ entry[signCol] in sigVal for entry in item if item[i][mapCol] in regions] )
        proteinsList = unique([ entry[proteinCol] for entry in in_region ])
        out.append([s_sign, s_sign_all, proteinsList])
    return out    

def getSignatureCountsFromDict(population, region, regions, mapCol, signCol, sigVal, baseCol, proteinCol):
    '''
    get signature counts from list of dicts
    return a list of lists indicating the number of signature mutations and number of -ve ctrl mutations, and the unique list of Uniprot IDs mutated
    '''
    if '_' in region:
        region = re.sub('_', '', region)
        regions = [re.sub('_', '', r) for r in regions]
    out = []
    for item in population:        
        s_sign = item[region + '_sig']
        s_sign_all = sum([item[r + '_sig'] for r in regions])
        proteinsList = unique(item['proteins'])
        out.append([s_sign, s_sign_all, proteinsList])
    return out    

def twoWayPvalFromNull(null_pop, val):
    percentile = scipy.stats.percentileofscore(null_pop, val)
    if percentile > 50.0:
        return 2 * ( 100 - percentile ) / 100.0
    else:
        return percentile * 2 / 100.0

def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]

def tidyProteinsList(proteinsList, maf):
    '''
    tidy uniprot IDs from MAF file by trying to use other mappings from other columns
    '''
    for i in range(len(proteinsList)):
        if proteinsList[i] == '':
            uniprot = maf[ i + 1 ][[entry for entry in range(len(maf[0])) \
                if maf[0][entry] == 'i_HGNC_UniProtIDsuppliedbyUniProt'][0]]
            if uniprot == '':
                uniprot = maf[ i + 1 ][[entry for entry in range(len(maf[0])) \
                    if maf[0][entry] == 'i_dbNSFP_Uniprot_acc'][0]]
                uniprot = uniprot.split(';')
                if len(uniprot) > 1:
                    uniprot_names = maf[ i + 1 ][[entry for entry in range(len(maf[0])) \
                        if maf[0][entry] == 'i_dbNSFP_Uniprot_acc'][0] + 1]
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

def tidySample(maf, countsTb, category, filter = False, geneColKey = 'Hugo_Symbol',\
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
    if category == 'structure':
        regions = ['surface', 'core', 'interact']
        mapType = 'struct_cat'
    if category == 'disorder':
        regions = ['intra_ord', 'intra_dis', 'inter_ord', 'inter_dis']
        mapType = 'disorder_cat'
        
    geneCol = [entry for entry in range(len(maf[0])) if maf[0][entry] == geneColKey][0]
    proteinCol = [entry for entry in range(len(maf[0])) if maf[0][entry] == proteinColKey][0]
    positionCol = [entry for entry in range(len(maf[0])) if maf[0][entry] == positionColKey][0]
    mapMafcol = [i for i in range(len(maf[0])) if maf[0][i] == mapType][0]

    if filter is not False:
        mafNew = [maf[0]]
        # take gene list in 'filter' and include only mutations which overlap with the list
        if 'NEGATIVE' in filter:
            for i in range(1, len(maf)):
                if maf[i][geneCol] not in filter:
                    mafNew.append(maf[i])
        else:
            for i in range(1, len(maf)):
                if maf[i][geneCol] in filter:
                    mafNew.append(maf[i])
        maf = mafNew

    proteinsList = [entry[proteinCol] for entry in maf[1:]]
    proteinsList = tidyProteinsList(proteinsList, maf)
        
    #remove those entries whose proteins/positions cannot map to zoomvar / structural regions
    for i in range(len(proteinsList)-1, -1, -1):
        maf[i+1][positionCol] = re.sub(" ", "", maf[i+1][positionCol])
        if maf[i+1][mapMafcol] == 'None':
            del proteinsList[i]
            del maf[i+1]
        else:
            maf[i+1][proteinCol] = proteinsList[i]
    proteinsList = unique(proteinsList)
    return [maf, proteinsList]        
    
def calculateEnrichment(maf, proteinsList, countsTb, category, sigVal, calc_bootstrap = True, calc_null = True, \
    proteinColKey = 'SwissProt_acc_Id', positionColKey = 'UniProt_AApos',\
    baseKey = 'Reference_Allele', altKey = "Tumor_Seq_Allele2", sigKey = 'APOBEC_mutation'):
    '''
    pipeline to generate:
    1) enrichment statistic for each region (type as defined by 'category')
    bootstrap enrichment statistic 95% CI
    2) p-value for the enrichment statistic derived from null distribution
    
    INPUT: tidied MAF
    '''
    if category == 'structure':
        regions = ['surface', 'core', 'interact']
        mapType = 'struct_cat'
    if category == 'disorder':
        regions = ['intra_ord', 'intra_dis', 'inter_ord', 'inter_dis']
        mapType = 'disorder_cat'                    
    proteinCol = [i for i in range(len(maf[0])) if maf[0][i] == proteinColKey][0]
    positionCol = [i for i in range(len(maf[0])) if maf[0][i] == positionColKey][0]
    mapCol = [i for i in range(len(maf[0])) if maf[0][i] == mapType][0]
    signCol = [i for i in range(len(maf[0])) if maf[0][i] == sigKey][0]
    baseCol = [i for i in range(len(maf[0])) if maf[0][i] == baseKey][0]
    altCol = [i for i in range(len(maf[0])) if maf[0][i] == altKey][0]                    
    binomial = {}
    enrich = {}
    bootstrap_ci = {}
    null_pval = {}
    null_ci = {}
    maf_sign_allregions = [maf[i] for i in range(len(maf[1:])) if maf[i][signCol].strip() in sigVal and  maf[i][mapCol].strip()!= 'None']#in regions ]
    mutlist = set()
    nmut_sign_allregions = 0
    #nmut_sign_allregions = sum([True for i in range(len(maf[1:])) if maf[i][signCol].strip() == '1' and maf[i][mapCol] in regions])
    for i in range(len(maf[1:])):
        if maf[i][proteinCol] + '_' + maf[i][positionCol] not in mutlist:
            if maf[i][signCol].strip() in sigVal and 'None' not in maf[i][mapCol].strip():
                nmut_sign_allregions += 1
            mutlist.add(maf[i][proteinCol].strip() + '_' + maf[i][positionCol].strip())
    print nmut_sign_allregions
    if calc_null:
        #calculate null population and derive empirical P-value
        pos_counts = collections.Counter([entry[proteinCol] for entry in maf_sign_allregions[1:]])
        if category == 'structure':
            null = simNull_struct(1000, pos_counts, countsTb)
        if category == 'disorder':
            null = simNull_disorder(1000, pos_counts, countsTb)
    if calc_bootstrap:
        #calculate bootstrap population and derive empirical CI
        bs = bootstrap(1000, nmut_sign_allregions, maf)
    for region in regions:
        nmut_sign = 0
        mutlist = set()
        #nmut_sign = sum([True for i in range(len(maf[1:])) \
        #                 if region in maf[i][mapCol] and maf[i][signCol].strip() in sigVal])
        for i in range(len(maf[1:])):
            if maf[i][proteinCol] + '_' + maf[i][positionCol] not in mutlist:
                if region in maf[i][mapCol].strip() and maf[i][signCol].strip() in sigVal:
                    nmut_sign += 1
                mutlist.add(maf[i][proteinCol].strip() + '_' + maf[i][positionCol].strip())
        print region, str(nmut_sign)
        #nmut_base = sum([True for i in range(len(maf[1:])) \
        #                if maf[i][mapCol] == region and \
        #                  ((maf[i][baseCol] == 'C' and maf[i][altCol] == 'T') or \
        #                   (maf[i][baseCol] == 'G' and maf[i][altCol] == 'A'))])
        proteins = unique([entry[proteinCol] for entry in maf_sign_allregions if entry[mapCol].strip() == region])
        #first calculate the observed statistic
        binomial[region] = binomCDF(proteins, countsTb, region, regions, \
                                    'sig', nmut_sign_allregions, nmut_sign)
        enrich[region] = apobecStructuralEnrichment(proteins, countsTb, region, regions, \
                                    nmut_sign, nmut_sign_allregions, context = 'sig')
        if calc_bootstrap:
        #calculate bootstrap population and derive empirical CI
            if nmut_sign > 0:
                # and nmut_base > 0:
                bootstrap_sign = getSignatureCountsFromMaf(bs, region, regions, mapCol, signCol, sigVal, baseCol, altCol, proteinCol)
                bootstrap_enrich = []
                for i in range(len(bootstrap_sign)):
                    bs_enrich = apobecStructuralEnrichment(bootstrap_sign[i][2], countsTb, region, regions, \
                        bootstrap_sign[i][0], bootstrap_sign[i][1], context = 'sig')
                    bootstrap_enrich.append(bs_enrich)
                bootstrap_enrich = [ entry for entry in bootstrap_enrich if entry is not None ]
                if len(bootstrap_enrich) > 0:
                    bootstrap_ci[region] = np.percentile(bootstrap_enrich, [2.5, 97.5])
                else:
                    bootstrap_ci[region] = ['None', 'None']
            else:
                bootstrap_ci[region] = ['None', 'None']
        else:
            bootstrap_ci[region] = ['None', 'None']

        if calc_null:
            null_sign = getSignatureCountsFromDict(null, region, regions, \
                                                   mapCol, signCol, sigVal, baseCol, proteinCol)
            null_enrich = [ apobecStructuralEnrichment(item[2], countsTb, region, regions, \
                item[0], item[1], context= 'sig') for item in null_sign]
            null_enrich = [entry for entry in null_enrich if entry is not None]
            if len(null_enrich) > 0:
                null_pval[region] = twoWayPvalFromNull(null_enrich, enrich[region])
                null_ci[region] = np.percentile(null_enrich, [2.5, 97.5])
            else:
                null_pval[region] = 'None'
                null_ci[region] = ['None', 'None']
        else:
            null_pval[region] = 'None'
            null_ci[region] = ['None', 'None']
    return {'binom_cdf': binomial, 'enrich': enrich, 'bootstrap_ci': bootstrap_ci, 'null_pval': null_pval, 'null_ci': null_ci}    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Script for computing the enrichment of mutations per sample and per cohort in each defined protein structural regions. Written for use with TCGA MAFs from Broad GDSC")
    parser.add_argument('--maf', help = 'mutation annotation file')
    parser.add_argument('--sigCol', help = 'column holding signature information')
    parser.add_argument('--sigVal', help = 'comma-separated values in sigCol corresponding to the signature of interest')
    parser.add_argument('--countsTb', help = 'per-protein counts of signature mutations at different regions (.csv)')
    parser.add_argument('--jobName', help = 'job name')
    parser.add_argument('--cohortfile', help = 'output filename (full path) for cohort-wide calculation')
    parser.add_argument('--caseList', help = 'list of cases to consider', default = None, required=False)
    parser.add_argument('--filter', help = 'whether to restrict analysis to a set of genes. Possible values = ["CGC", "BaileyEtAl_pancancer", "BaileyEtAl_ctSpecific", "all"]', default = 'all')
    args = parser.parse_args()

    print args.filter
    if args.filter == 'CGC':
        # read the list of genes to analyse
        drivers = set() # for all drivers
        with open('CGC_COSMICv86_20181030_Tier1.tsv', 'r') as instream:
            read = csv.reader(instream, delimiter="\t")
            for row in read:
                drivers.add(row[0])
        drivers = list(drivers)
    elif args.filter == 'nonCGC':
        # read the list of genes to analyse
        drivers = set() 
        drivers.add('NEGATIVE')
        with open('CGC_COSMICv86_20181030_Tier1.tsv', 'r') as instream:
            read = csv.reader(instream, delimiter="\t")
            for row in read:
                drivers.add(row[0])
        drivers = list(drivers)
    elif args.filter == 'BaileyEtAl_pancancer':
        drivers = set()
        with open('BaileyEtAl_CancerDrivers.txt', 'r') as instream:
            read = csv.reader(instream, delimiter="\t")
            for row in read:
                if row[1] == 'PANCAN':
                    drivers.add(row[0])     # for all drivers
        drivers = list(drivers)
    elif args.filter == 'BaileyEtAl_ctSpecific':
        drivers = list()
        with open('BaileyEtAl_CancerDrivers.txt', 'r') as instream:
            read = csv.reader(instream, delimiter="\t")
            for row in read:
                if row[1] in ['PANCAN', args.jobName]: # for cancer specific drivers
                    drivers.append(row[0]) # for cancer specific drivers
    elif args.filter == 'all':
        drivers = False

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
    print "Analysing cohort: " + args.jobName        
    sigs = args.sigVal.split(',')

    # get only the indicated cases:
    if args.caseList is not None:
        with open(args.caseList, 'r') as instream:
            filterSamples = instream.readlines()
            filterSamples = [l[:-1] for l in filterSamples]
        print len(filterSamples)
    sampleCol = [entry for entry in range(len(maf[0])) \
                 if maf[0][entry] == 'Tumor_Sample_Barcode'][0]
    if args.caseList is not None:
        mafSelected = [entry for entry in maf if entry[sampleCol] == "Tumor_Sample_Barcode" or \
            entry[sampleCol][:12] in filterSamples]
        maf = mafSelected
    samples = list(set([ entry[sampleCol][:16] for entry in maf[1:] ]))
    print(len(samples))
    variantCol = [entry for entry in range(len(maf[0])) \
                  if maf[0][entry] == 'Variant_Classification'][0]
    wtCol = [entry for entry in range(len(maf[0])) if maf[0][entry] == 'origAA'][0]
    mutCol = [entry for entry in range(len(maf[0])) if maf[0][entry] == 'mutAA'][0]
    maf = [maf[0]] + [entry for entry in maf if entry[wtCol] in amino_acids and \
                      entry[mutCol] in amino_acids and entry[wtCol] != entry[mutCol] ]
    #tidy the MAF
    print "Tidying MAF file ..."
    maf_struct = tidySample(maf, countsTb, 'structure', filter = drivers)
    print len(maf_struct[0])
    #cohort-wide analysis first
    with open(args.cohortfile, 'a') as outstream:
        wr = csv.writer(outstream)        
        print "Now calculate enrichment ..."
        print "Structural regions ..."
        out_struct = calculateEnrichment(maf_struct[0], maf_struct[1], \
            countsTb, 'structure', sigVal = sigs, sigKey = args.sigCol, calc_bootstrap = True, calc_null = True)
        wr.writerow([args.jobName, \
            len(maf_struct[0]) - 1, \
            str(out_struct['binom_cdf']['surface']), str(out_struct['enrich']['surface']), \
            out_struct['bootstrap_ci']['surface'][0], \
            out_struct['bootstrap_ci']['surface'][1], \
            out_struct['null_ci']['surface'][0], \
            out_struct['null_ci']['surface'][1], \
            out_struct['null_pval']['surface'], \
            str(out_struct['binom_cdf']['core']), str(out_struct['enrich']['core']), \
            out_struct['bootstrap_ci']['core'][0], \
            out_struct['bootstrap_ci']['core'][1], \
            out_struct['null_ci']['core'][0], \
            out_struct['null_ci']['core'][1], \
            out_struct['null_pval']['core'], \
            str(out_struct['binom_cdf']['interact']), \
            str(out_struct['enrich']['interact']), \
            out_struct['bootstrap_ci']['interact'][0], \
            out_struct['bootstrap_ci']['interact'][1], \
            out_struct['null_ci']['interact'][0], \
            out_struct['null_ci']['interact'][1], \
            out_struct['null_pval']['interact']])

