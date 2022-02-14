import pandas as pd
import numpy as np
import re

# functions for one-hot encoding

def decomposeMutSig(mutsig):
    """
    given a mutsig string of format 'W[X>Y]Z', decompose into position -1, +1 and the actual substitution
    Input:
        mutsig: string of format 'W[X>Y]Z'
    Output:
        tuple of (pos-1, pos+1, mut). So 'W[X>Y]Z' --> ('W', 'Z', 'X>Y')
    """
    pos_minus1 = mutsig[0]
    pos_plus1 = mutsig[-1]
    mut = re.search('[AGCT]>[AGCT]', mutsig).group(0)
    return (pos_minus1, pos_plus1, mut)
    

def onehot(data, levels):
    """
    one-hot encoding of the 'data'.
    Input:
        data - the list/pd.Series of data to be encoded
        levels - list of possible 'factor' levels of this list/series
    Output:
        np.array of shape (n_data, n_levels)
    """
    encode = np.zeros((len(data), len(levels)))
    for i, val in enumerate(data):
        encode[i, levels.index(val)] = 1
    return encode    

def prepare_data(file, delimiter, output_col, WT_AA = True, MUT_AA = True, \
                Pos_1 = True, Subs = True, PhyloP = True, QSASA = True):
    """
    wrapper function to read in the data and do required preprocessing
    Input:
        file - str, filepath. Has header in the first column and no row index columns.
        delimiter - str, delimiter used in 'file'
        output_col, str, column name in `file` corresponding to the output labels 'Y' 
    Output: 
        tuple of 2 or 3 elements:
        - np.array of shape (db.shape[0], 58) where db is the resulting data frame from reading in `file`. These are the features 'X'.
        - np.array of shape (db.shape[0], 1). These are the output 'Y' to be predicted.
        - np.array of shape (db.shape[0], 1). **Return only if 'scoreset' is present in the columns of db.** The labels of the scoreset of each obs. For sampling to ensure equal ratio of labels individually in each score set.
    """
    db = pd.read_csv(file, sep = delimiter, header=0, index_col=False, na_values = 'None')
    print('Data read in. Now process them')
    if PhyloP:
        # normalisation: PhyloP score columns / 10 so that it ranges from -1 to 1.
        db[['PhyloP_-1', 'PhyloP_0', 'PhyloP_1']] /= 10
        print('PhyloP score normalised to [-1, 1].')

    # decomposeMutSig on 'MutSig'
    db[['Pos_-1', 'Pos_+1', 'Mut']] = pd.DataFrame([decomposeMutSig(val) for val in db['MutSig']], index = db.index)
    print('MutSig parsed.')

    # one-hot on the new columns, WT_AA, MUT_AA.
    AAs = list('ACDEFGHIKLMNPQRSTVWY')
    feature_cols = []
    if WT_AA:
        print('WT_AA one-hot encoding ...')
        db[ ['WT_AA_' + aa for aa in AAs] ] = pd.DataFrame(onehot(db['WT_AA'], AAs), index = db.index)
        feature_cols +=  ['WT_AA_' + aa for aa in AAs]
    if MUT_AA:
        print('MUT_AA one-hot encoding ...')
        db[ ['MUT_AA_' + aa for aa in AAs] ] = pd.DataFrame(onehot(db['MUT_AA'], AAs), index = db.index)
        feature_cols +=  ['MUT_AA_' + aa for aa in AAs]
    bs = list('AGCT')
    if Pos_1:
        print('-1 and +1 DNA pos one-hot encoding ...')
        db[ ['Pos_-1_' + b for b in bs] ] = pd.DataFrame(onehot(db['Pos_-1'], bs), index = db.index)
        db[ ['Pos_+1_' + b for b in bs] ] = pd.DataFrame(onehot(db['Pos_+1'], bs), index = db.index)
        feature_cols +=  ['Pos_-1_' + b for b in bs]
        feature_cols +=  ['Pos_+1_' + b for b in bs]
    subs = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    if Subs:    
        print('DNA substitution one-hot encoding ...')
        db[ ['Mut_' + m for m in subs] ] = pd.DataFrame(onehot(db['Mut'], subs), index = db.index)
        feature_cols +=  ['Mut_' + m for m in subs] 
    if PhyloP:
        feature_cols += ['PhyloP_0']
        if Pos_1:
            feature_cols +=  ['PhyloP_-1', 'PhyloP_1']
    if QSASA:
        db['Q(SASA)'] -= 0.15
        feature_cols += ['Q(SASA)']
    if 'scoreset' in db.columns:
        return np.array( db[ feature_cols] ), np.array( db[output_col] ).reshape((-1,1)), np.array(db['scoreset']).reshape((-1, 1))
    else:
        return np.array( db[ feature_cols] ), np.array( db[output_col] ).reshape((-1,1))
