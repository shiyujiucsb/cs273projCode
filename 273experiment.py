#!/usr/bin/python3
#
# Authors: Shiyu Ji and Kun Wan
# A new bound on worst case error for FI mining approximation.
# 2016/10/17
#

'''
Exact solutions.
Input: dataset file name, # itemsets
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
Input: dataset file name, # itemsets, # samples.
Output: list of approximated itemset frequencies
'''
def apprxFreqs(fname, nItemsets, nSamples):
    freqs = [0 for _ in range(nItemsets+1)]
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    nTransactions = len(lines)
    from random import randint
    for _ in range(nSamples):
        line = lines[randint(0,nTransactions-1)]
        transaction = [int(i) for i in line.strip().split(' ')]
        for itemset in transaction:
            if 0 <= itemset <= nItemsets:
                freqs[itemset] += 1
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
    
'''
Approximation with new bound.
Input: dataset file name, # itemsets, epsilon, delta
Output: approximated frequencies, # samples
'''
def newBoundApprox(fname, nItemsets, epsilon, delta):
    f = open(fname, 'r')
    nTransactions = 0
    for _ in f: nTransactions += 1
    f.close()
    
    from math import log
    n = int((log(2*nItemsets)-log(delta))/(2*epsilon*epsilon))+1
    if n > nTransactions:
        return calcExactFreqs(fname, nItemsets), nTransactions
    return apprxFreqs(fname, nItemsets, n), n
    
'''
Classical approximation using old bound (RU15@KDD)
Input: dataset file name, # itemsets, epsilon, delta
Output: approximations, # samples
'''
def RUApprox(fname, nItemsets, epsilon, delta):
    freqs = [0 for _ in range(nItemsets+1)]
    f = open(fname, 'r')
    lines = f.readlines()
    nTransactions = len(lines)
    f.close()
    from random import randint
    from math import log
    n = 0
    dn = int(2*log(2/delta)/(epsilon*epsilon))+1
    while True:
        for _ in range(n, n+dn):
            line = lines[randint(0,nTransactions-1)]
            transaction = [int(i) for i in line.strip().split(' ')]
            for itemset in transaction:
                if 0 <= itemset <= nItemsets:
                    freqs[itemset] += 1
        n = n+dn
        ell = max(freqs)**0.5
        Delta = 2*ell/n*(2*log(nItemsets))**0.5 + (2*log(2/delta)/n)**0.5
        if Delta <= epsilon or n >= nTransactions:
            break
        dn = int(2*log(2/delta)/((epsilon-2*ell/n*(2*log(nItemsets))**0.5)**2))+1-n
        dn = max(dn,1)
    return list(map(lambda x:x*1.0/n,freqs)), n
