from collections import Counter
import numpy as np

def getKmers(string,kV):
    return [ string[i:i+kV] for i in range(0,len(string)-kV+1)]

def getCounters(seqs,kmer,goodIdx):
    working = seqs[goodIdx]
    # transform the dna into their kmers. Use a counter structure, 
    # which is much faster. it stored the unique elements and counts
# see: https://docs.python.org/2/library/collections.html#collections.Counter
    return map(lambda x: Counter(getKmers(x,kmer)),working)

def getMinKFromCounters(counters,kmer,goodIdx,returnKmers=None,
                        printProgress=True):
    # assuming 'counters' is a counter object with each {kmer : count} as a dict
    kmerDat = [ tmpSet.most_common() for tmpSet in counters ]
    kmerSet = map(lambda x: [item[0] for item in x],kmerDat)
    kmerCounts = map(lambda x: [item[1] for item in x],kmerDat)
    # maximum count (ie: the maximum number of times *any* kmer happens
    # this needs to be 1 to find the appropriate value...
    maxCount = map(max,kmerCounts)
    # get the index where the max was one (ie: at most 1 of the kmers)
    bestIdx = [ goodIdx[i] for i,c in enumerate(maxCount)  if c<=1]
    # update the bookkeeping stuff; of the originals (goodIdx.size), how many
    # did we find the right kmer for len(bestIdx)
    kNotFoundNum = goodIdx.size - len(bestIdx)
    if (printProgress):
        # assume we are just interested in the first kmer set, for debugging..
        print("GetMinK: {:d} seqs left, average max of {:.1f} {:d}-mers".
              format(kNotFoundNum,np.mean(maxCount),kmer))
    kmer += 1
    if (returnKmers is None):
        return kmer,kNotFoundNum,bestIdx
    else:
        # also return the data the the 'returnKmers' most common ones.
        return kmer,kNotFoundNum,bestIdx,kmerDat[0][:returnKmers]

def getMinKIdxAndCount(seqs,kmer,goodIdx,returnKmers=None,printProgress=True):
    # returnKmers: if not none, returns this many (including -1, all)
    # of the kmers found 
    # get the set of each of the kmers, using the unique set of the counters
# see: https://docs.python.org/2/library/collections.html#collections.Counter
    counters = getCounters(seqs,kmer,goodIdx)
    # For each counter, get all its elements (most common defaults to n)
    # and record the count (second element, first index)
# see: https://docs.python.org/2/library/collections.html#collections.Counter
    return getMinKFromCounters(counters,kmer,goodIdx,returnKmers,printProgress)
