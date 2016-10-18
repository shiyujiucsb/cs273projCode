#!/usr/bin/python3
#
# Authors: Shiyu Ji and Kun Wan
# A new bound on worst case error for FI mining approximation.
# 2016/10/17
#

'''
Exact solutions.
Input: data set file name, # itemsets
Output: list of itemset frequencies
'''
def calcExactFreqs(fname, nItemsets):
    freqs = [0 for _ in range(nItemsets+1)]
    f = open(fname, 'r')
    nTransactions=0
    for line in f:
        nTransactions+=1
        transaction = [int(i) for i in line.strip().split(' ')]
        for itemset in transaction:
            if 0 <= itemset <= nItemsets:
                freqs[itemset] += 1
    f.close()
    return list(map(lambda x:x*1.0/nTransactions, freqs))

'''
Approximated solutions given sample size.
Input: data set file name, # itemsets, # samples.
Output: list of approximated itemset frequencies
'''
def apprxFreqs(fname, nItemsets, nSamples):
    freqs = [0 for _ in range(nItemsets+1)]
    f = open(fname, 'r')
    lines = f.readlines()
    nTransactions = len(lines)
    from random import randint
    for _ in range(nSamples):
        line = lines[randint(1,nTransactions)]
        transaction = [int(i) for i in line.strip().split(' ')]
        for itemset in transaction:
            if 0 <= itemset <= nItemsets:
                freqs[itemset] += 1
    f.close()
    return list(map(lambda x:x*1.0/nSamples, freqs))

'''
Find the maximal error (worst case).
Input: two freq arrays with identical lengths
Output: the maximal error
'''
def worstError(freq1, freq2):
    if len(freq1) != len(freq2) or not freq1 or not freq2:
        return None
    return max([abs(freq1[i]-freq2[i]) for i in range(len(freq1))])
    
