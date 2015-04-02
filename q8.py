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

def getRandomDNA(chars,dnaLen,numToChoose,weights):
    rawArr = np.random.choice(chars,size=(numToChoose,dnaLen))
    # POST: rawArr[i][j] has character [j] of dna string [i] 
    # (ie: each row is a string)
    # combine all the columns in a given row for the actual strings
    return np.array([ "".join(r) for r in rawArr],dtype=np.object)
def getKmers(string,k):
    return [ string[i:i+k] for i in range(0,len(string)-k+1)]
    
def getMinK(seqs):
    arrK =np.ones(numOligos,dtype=np.uint32) * -1
    kNotFoundNum = arrK.size
    goodIdx =np.where(arrK < 0)[0]
    kmer = 1
    while (kNotFoundNum > 0):
        working = seqs[goodIdx]
        # transform the dna into their kmers
        kmers = map(lambda x: getKmers(x,kmer),working)
        # get the count of each kmers occurance
        kmerSet = map(set,kmers)
        # For each set of kmers kmerSet[i], count occurences in kmers[i]
        kmerCounts = [ [ kmers[i].count(c) for c in tmpSet ] \
                       for i,tmpSet in enumerate(kmerSet) ]
        # get the maximum count 
        maxCount = map(max,kmerCounts)
        # get the index where the max was one (ie: at most 1 of the kmers)
        bestIdx = [ goodIdx[i] for i,c in enumerate(maxCount)  if c==1]
        arrK[bestIdx] = kmer
        goodIdx =np.where(arrK < 0)[0]
        kNotFoundNum = goodIdx.size
        print(kNotFoundNum)
        kmer += 1
    return arrK

def getKSequence(lenArr,numOligos,weights,chars):
    toRet = []
    for idx,lenV in enumerate(lenArr):
        dna = getRandomDNA(chars,lenV,numOligos,weights)
        toRet.append(getMinK(dna))
        print("{:d}/{:d}".format(idx+1,len(lenArr)))
    return toRet

def plotAll(kArrs,outDir):
    maxKeach = [max(k) for k in kArrs ]
    maxK = max(maxKeach)
    bins = range(maxK+1)
    numTrials = len(kArr)
    means  = np.array([np.mean(k) for k in kArrs])
    stdevs = np.array([np.std(k) for k in kArrs])
    for i,k in enumerate(kArrs):
        fig = pPlotUtil.figure()
        plt.hist(k,bins=bins,align='left',label='Data from {:d} sequences'.
                 format(int(numOligos)))
        mean = means[i]
        plt.axvline(mean,color='r',label="Mean:{:.3f}".format(mean),
                    linewidth=2.0)
        plt.xlim([0,maxK])
        plt.xlabel('K, minimum k-mer with at most 1 occurence in DNA')
        plt.ylabel('Number of occurences')
        plt.title('K histogram for a dna of length {:d}'.format(lengths[i]))
        plt.legend()
        pPlotUtil.savefig(fig,outDir + "k{:d}".format(i))
    return means,stdevs

table = [ ['A','G','C','T'],
          [0.1,0.4,0.2,0.3]]
chars = table[0]
weights = table[1]
q = max(weights)
 # use for easy python string generation
# POST: bag has the characters in their appropriate proportion
numOligos = 10e3
lengths = np.array([32,64,128,256,512])
# save the K array: minimum k to have at most one k-mer
# initialize to -1, so that we know when we have the minimum
outDir = "./out/"
pGenUtil.ensureDirExists(outDir)

kArr = pCheckUtil.getCheckpoint('./tmp/check.pkl',getKSequence,False,
                                lengths,numOligos,weights,chars)
meanVals,std = pCheckUtil.getCheckpoint('./tmp/meanStd.pkl',plotAll,False,kArr,
                                    outDir)

fig = pPlotUtil.figure()
plt.errorbar(x=lengths,y=meanVals,yerr=std,fmt='ro-',label='Mean K')
plt.plot(lengths,np.log(lengths+1)/np.log(1/q),'b--',label='Log_[1/q](l+1)')
plt.xlabel('DNA Length (l)')
plt.ylabel('Mean K value')
plt.title('Mean K vs length')
plt.legend(loc='best')
pPlotUtil.savefig(fig,outDir + 'k_v_len')
