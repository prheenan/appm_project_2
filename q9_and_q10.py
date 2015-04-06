# force floating point division. Can still use integer with //
from __future__ import division
# This file is used for importing the common utilities classes.
import numpy as np
import matplotlib.pyplot as plt
# need to add the utilities class. Want 'home' to be platform independent
from os.path import expanduser
home = expanduser("~")
# get the utilties directory (assume it lives in ~/utilities/python)
# but simple to change
path= home +"/utilities/python"
import sys
sys.path.append(path)
# import the patrick-specific utilities
import GenUtilities  as pGenUtil
import PlotUtilities as pPlotUtil
import CheckpointUtilities as pCheckUtil
from util import getMinKIdxAndCount

baseDict = {'A':'T',
            'T':'A',
            'C':'G',
            'G':'C'}

def getRevComplement(seq):
    return ''.join([ baseDict[base] for base in seq][::-1])

def getProportions(seq,revComp):
    # just concatenate the sequence and its rev comp for counting purposes
    psuedoSeq = seq + revComp
    bases = baseDict.keys()
    nBases = len(psuedoSeq)
    # base proportions array: count the number of times a base occurs, divide
    # by the total number of bases
    baseProps = map(lambda x: psuedoSeq.count(x)/nBases,bases)
    return np.array(baseProps,dtype=np.float64)

def printProportions(baseProps):
    for proportion,base in zip(baseProps,baseDict.keys()):
        print("The proportion of {:s} is {:.4f}".format(base,proportion))
    print('\tSums to... {:.5f} == 1 (I hope!)'.format(sum(baseProps)))

def findKmers(sequence,seqComplement):
    kmer = 1
    circular = lambda seq ,k: seq + seq[:(k-1)] if k > 1 else seq
    goodIdx =np.array([0])
    kNotFound = 1
    avgKmerCountS1 = []
    avgKmerCountS2 = []
    # kmers to save each round
    topKmers = 10
    # maximum numbers of rounds to save
    maxK = 300
    # for each round [i] and both sequences [j],  store the top kmers [k]
    saveKs = np.empty((maxK,2,topKmers),dtype=np.object)
    while (kNotFound >0):
        s1 = circular(sequence,kmer)
        s2 = circular(seqComplement,kmer)
        # need to pass an array to minK..
        xx,notFoundS1,xx,commonKmers1 = getMinKIdxAndCount([[s1]],kmer,goodIdx,
                                                          topKmers)
        xx,notFoundS2,xx,commonKmers2 = getMinKIdxAndCount([[s2]],kmer,goodIdx,
                                                          topKmers)
        saveKs[kmer-1,0,:min(topKmers,len(commonKmers1))] = commonKmers1
        saveKs[kmer-1,1,:min(topKmers,len(commonKmers2))] = commonKmers2
        # continue going until both the forward and reverse sequences satistfy
        kNotFound = max(notFoundS1,notFoundS2)
        # keep arrK at -1 if either sequences don't have the K yet.
        kmer += 1
    # go up to but not including the kmer we found; this gets up the last 
    # k for which we had more than one k-mer in the sequence
    return kmer,saveKs[:kmer-1,:,:]

if __name__ == '__main__':

    sequenceFile = "./data/GCA_000027325.1_ASM2732v1_genomic_sequence.txt"
    with open(sequenceFile) as f:
        sequence = f.read().strip().upper()
    # get the reverse complement
    seqComplement = getRevComplement(sequence)
    baseProps = pCheckUtil.getCheckpoint('./tmp/proportions.npz',\
                            getProportions,False,sequence,seqComplement)
    printProportions(baseProps)
    # lazy way of making the sequence circularized
    kmer,saveKs = pCheckUtil.getCheckpoint('./tmp/minK.pkl',findKmers,
                                           False,sequence,seqComplement)

