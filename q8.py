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

def getRandomDNA(chars,dnaLen,numToChoose,weights):
    rawArr = np.random.choice(chars,size=(numToChoose,dnaLen),p=weights)
    # POST: rawArr[i][j] has character [j] of dna string [i] 
    # (ie: each row is a string)
    # combine all the columns in a given row for the actual strings
    return np.array([ "".join(r) for r in rawArr],dtype=np.object)

def getMinK(seqs,printProgress=False):
    # get the minimum integer k such that each -mer occurs once in each of seqs
    # note: we will use '-1' to note 'we haven't found the minimum k yet'
    numOligos = len(seqs)
    arrK =np.ones(numOligos,dtype=np.uint32) * -1
    # use the following two bookkeeping variables to help with looping
    kNotFoundNum = arrK.size
    goodIdx =np.where(arrK < 0)[0]
    kmer = 1
    while (kNotFoundNum > 0):
        kmer,kNotFoundNum,bestIdx = getMinKIdxAndCount(seqs,kmer,goodIdx,
                                    printProgress=printProgress)
        arrK[bestIdx] = kmer
        goodIdx =np.where(arrK < 0)[0]
    return arrK

def getKSequence(lenArr,numOligos,weights,chars):
    toRet = []
    for idx,lenV in enumerate(lenArr):
        # get the random DNA according to the lengths we have, the chars/weights
        dna = getRandomDNA(chars,lenV,numOligos,weights)
        # append the minimum k such that each kmer appears at most once
        toRet.append(getMinK(dna))
        # progress bar! 
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
                 format(int(numOligos)),normed=True)
        mean = means[i]
        plt.axvline(mean,color='r',label="Mean:{:.3f}".format(mean),
                    linewidth=2.0)
        plt.xlim([0,maxK])
        plt.xlabel('K, minimum k-mer with at most 1 occurence in DNA sequence')
        plt.ylabel('Proportion of occurences')
        plt.title('K histogram (normalized) for DNA sequences  of length {:d}'.
                  format(lengths[i]))
        plt.legend()
        pPlotUtil.savefig(fig,outDir + "k{:d}".format(i))
    return means,stdevs

def getTheoryK(lengths,q):
    c0 = np.log(1/q)
    return np.log(lengths + 1)/c0


def testDnaGeneration(chars,lengths,numOligos,weights):
    numChars = len(chars)
    oligoWeights = np.zeros((len(lengths),numChars,numOligos))
    for i,l in enumerate(lengths):
        dna = getRandomDNA(chars,l,numOligos,weights)
        for j,c in enumerate(chars):
            oligoWeights[i][j][:] = [ d.count(c)/l for d in dna ]
    for k,c in enumerate(chars):
        print("testDNA, expected/actual %{:s} is {:.3f}/{:.3f}".
              format(c,weights[k],np.mean(oligoWeights[:,k,:])))

def plotError(expected,actual,xV,xlab,ylab,title,ax,relative=False):
    plt.title(title)
    delta = np.abs(expected-actual)
    y = delta/np.minimum(expected,actual) if relative else delta
    plt.plot(xV,y,'r-')
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    ax.set_xscale('log')

if __name__ == '__main__':
    table = [ ['A','G','C','T'],
              [0.1,0.4,0.2,0.3]]
    chars = table[0]
    weights = table[1]
    q = max(weights)
    # use for easy python string generation
    numOligos = 10e3
    lengths = np.array([2,4,8,16,32,64,128,256,512])
    # save the K array: minimum k to have at most one k-mer
    # initialize to -1, so that we know when we have the minimum
    outDir = "./out/"
    pGenUtil.ensureDirExists(outDir)
    forceRun = False
    test = False
    # use checkpointing to save data, since it takes forever
    kArr = pCheckUtil.getCheckpoint('./tmp/check.pkl',getKSequence,forceRun,
                                    lengths,numOligos,weights,chars)
    meanVals,std = pCheckUtil.getCheckpoint('./tmp/meanStd.pkl',plotAll,
                                            forceRun,kArr,outDir)
    if (test):
        testDnaGeneration(chars,lengths,numOligos,weights)
    # plot the mean k vs dna length, l (in theory, k is approx log_1/q(l+1))
    fig = pPlotUtil.figure()
    ax = plt.subplot(1,3,1)
    plt.errorbar(x=lengths,y=meanVals,yerr=std,fmt='ro-',label='Mean K')
    tKVals = getTheoryK(lengths,q)
    plt.plot(lengths,tKVals,'b--',label='Log_[1/q](l+1)')
    xLab = 'DNA Length (l)'
    plt.xlabel(xLab)
    plt.ylabel('Mean K value')
    plt.title('Mean K vs length')
    ax.set_xscale('log')
    plt.legend(loc='best')
    ax = plt.subplot(1,3,2)
    plotError(meanVals,tKVals,lengths,xLab,'Absolute Error in Mean K ',
              'Absolute error in Mean K',ax,relative=False)
    ax = plt.subplot(1,3,3)
    plotError(meanVals,tKVals,lengths,xLab,'Relative Error in Mean K [0-->1]',
              'Relative error in Mean K',ax,relative=True)
    pPlotUtil.savefig(fig,outDir + 'k_v_len')

