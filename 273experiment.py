#!/usr/bin/python3
#
# Authors: Shiyu Ji and Kun Wan
# A new bound on worst case error for FI mining approximation.
# 2016/10/22
#

'''
Determine whether an itemset is in the transaction
Input: a transaction and an itemset
Output: True if the transaction contains the itemset
'''
def isContained(itemset, transaction):
    transactionSet = set(transaction)
    for item in itemset:
        if item not in transactionSet:
            return False
    return True

'''
Exact solutions.
Input: dataset file name, itemsets
Output: list of itemset frequencies
'''
def calcExactFreqs(fname, itemsets):
    freqs = [0 for _ in range(len(itemsets))]
    f = open(fname, 'r')
    nTransactions=0
    for line in f:
        nTransactions+=1
        transaction = [int(i) for i in line.strip().split(' ')]
        for i in range(len(itemsets)):
            if isContained(itemsets[i], transaction):
                freqs[i] += 1
    f.close()
    return list(map(lambda x:x*1.0/nTransactions, freqs))

'''
Approximated solutions given sample size.
Input: dataset file name, itemsets, # samples.
Output: list of approximated itemset frequencies
'''
def apprxFreqs(fname, itemsets, nSamples):
    freqs = [0 for _ in range(len(itemsets))]
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    nTransactions = len(lines)
    from random import randint
    for _ in range(nSamples):
        line = lines[randint(0,nTransactions-1)]
        transaction = [int(i) for i in line.strip().split(' ')]
        for i in range(len(itemsets)):
            if isContained(itemsets[i], transaction):
                freqs[i] += 1
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
Input: dataset file name, itemsets, epsilon, delta
Output: approximated frequencies, # samples
'''
def newBoundApprox(fname, itemsets, epsilon, delta):
    f = open(fname, 'r')
    nTransactions = 0
    for _ in f: nTransactions += 1
    f.close()
    
    from math import log
    n = int((log(2*len(itemsets))-log(delta))/(2*epsilon*epsilon))+1
    if n > nTransactions:
        return calcExactFreqs(fname, itemsets), nTransactions
    return apprxFreqs(fname, itemsets, n), n
    
'''
Classical approximation using old bound (RU15@KDD)
Input: dataset file name, itemsets, epsilon, delta
Output: approximations, # samples
'''
def RUApprox(fname, itemsets, epsilon, delta):
    freqs = [0 for _ in range(len(itemsets))]
    f = open(fname, 'r')
    lines = f.readlines()
    nTransactions = len(lines)
    f.close()
    from random import randint
    from math import log
    n = 0
    dn = min(int(2*log(2/delta)/(epsilon*epsilon))+1, nTransactions)
    while n < nTransactions:
        for _ in range(n, n+dn):
            line = lines[randint(0,nTransactions-1)]
            transaction = [int(i) for i in line.strip().split(' ')]
            for i in range(len(itemsets)):
                if isContained(itemsets[i], transaction):
                    freqs[i] += 1
        n = n+dn
        ell = max(freqs)**0.5
        Delta = 2*ell/n*(2*log(len(itemsets)))**0.5 + (2*log(2/delta)/n)**0.5
        if Delta <= epsilon or n >= nTransactions:
            break
        dn = int(2*log(2/delta)/((epsilon-2*ell/n*(2*log(len(itemsets)))**0.5)**2))+1-n
        dn = max(dn,1)
    return list(map(lambda x:x*1.0/n,freqs)), n

'''
Exact top-k algorithm.
Input: dataset file name, itemsets, top-K
Output: top-K itemsets and their frequencies
'''
def exactTopK(fname, itemsets, K):
    freqs = calcExactFreqs(fname, itemsets)
    topKFIs = sorted(range(len(itemsets)), key=lambda x:-freqs[x])
    return list(map(lambda x:str(itemsets[x]), topKFIs)), sorted(freqs, key=lambda x:-x)

'''
Approx top-k algorithm with new bound.
Input: dataset file name, itemsets, top-K, epsilon, delta
Output: top-K itemsets and their frequencies
'''
def newBoundTopK(fname, itemsets, K, epsilon, delta):
    freqs, n = newBoundApprox(fname, itemsets, epsilon, delta)
    topKFIs = sorted(range(len(itemsets)), key=lambda x:-freqs[x])
    return list(map(lambda x:str(itemsets[x]), topKFIs)), sorted(freqs, key=lambda x:-x)

'''
Approx top-k algorithm by RU's classical method.
Input: dataset file name, itemsets, top-K, epsilon, delta
Output: top-K itemsets and their frequencies
'''
def RUtopK(fname, itemsets, K, epsilon, delta):
    freqs, n = RUApprox(fname, itemsets, epsilon, delta)
    topKFIs = sorted(range(len(itemsets)), key=lambda x:-freqs[x])
    return list(map(lambda x:str(itemsets[x]), topKFIs)), sorted(freqs, key=lambda x:-x)

'''
Our top-K approximate algorithm
Input: dataset file name, itemsets, top-K, sample increase, delta
    max sample ratio alpha
Output: top-K itemsets and their frequencies, # samples
'''
def topK(fname, itemsets, K, sampleInc, delta, alpha):
    freqs = [0 for _ in range(len(itemsets))]
    f = open(fname, 'r')
    lines = f.readlines()
    nTransactions = len(lines)
    f.close()
    from random import randint
    from math import log, exp
    n = 0
    while True:
        for _ in range(sampleInc):
            line = lines[randint(0,nTransactions-1)]
            transaction = [int(i) for i in line.strip().split(' ')]
            for i in range(len(itemsets)):
                if isContained(itemsets[i], transaction):
                    freqs[i] = freqs[i]+1
        n += sampleInc
        tmp = sorted(freqs, key=lambda x:-x)
        middle = (tmp[K-1]+tmp[K])*.5
        errorProb=0.0
        for i in range(len(tmp)):
            dProb = exp(-2.0/n*abs(tmp[i]-middle)**2)
            errorProb +=dProb
            if i>=K:
                if errorProb>delta or errorProb+(len(tmp)-1-i)*dProb<=delta:
                    break
        if errorProb <= delta or n>=alpha*nTransactions:
            return list(map(lambda x:itemsets[x], \
                            sorted(range(len(itemsets)), key=lambda x:-freqs[x])[:K])),\
                    list(map(lambda x:x*1.0/n, tmp[:K])), n

'''
A-Priori Top-K algorithm
Input: dataset file name, # items, top-K, sample increase, delta,
    max sample ratio alpha
Output: top-K itemsets and their frequencies, # samples
'''
def AprioriTopK(fname, I, K, sampleInc, delta, alpha):
    if I*I<K:
        print('Too few items for K')
        return
    singletons = [[i] for i in range(1, I+1)]
    doubletons = []
    if K>=I:
        topKSingle = singletons
    else:
        topKSingle, topKSingleFreqs, nSmp = topK(fname, singletons, K, sampleInc, delta, alpha)
    topKSingle.sort()
    for i in range(min(K,I)):
        for j in range(i+1, min(K,I)):
            doubletons.append(topKSingle[i] + topKSingle[j])
    topKDouble, topKDoubleFreqs, nSmp = topK(fname, doubletons, K, sampleInc, delta, alpha)
    return topKDouble, topKDoubleFreqs, nSmp

'''
Precision test (epsilon-based): select K most frequent pairs. 
Input: dataset file name, # items, top-K, epsilon, delta
Output: precision, running time.
'''
def testEpsilonPrecision(fname, I, K, epsilon, delta):
    if I*I<K:
        print('Too few items for K')
        return
    singletons = [[i] for i in range(1, I+1)]
    doubletons = []
    if K>=I:
        topKSingle = singletons
    else:
        topKSingle, topKSingleFreqs = exactTopK(fname, singletons, K)
    topKSingle.sort()
    for i in range(min(K,I)):
        for j in range(i+1,min(K,I)):
            doubletons.append(topKSingle[i] + topKSingle[j])
            
    from time import time
    startTimeExact = time()
    exactFIs, exactFreqs = exactTopK(fname, doubletons, K)
    startTimeOurs = endTimeExact = time()
    approxFIs, approxFreqs = newBoundTopK(fname, doubletons, K, epsilon, delta)
    startTimeRU = endTimeOurs = time()
    RUapprxFIs, RUapprxFreqs = RUtopK(fname, doubletons, K, epsilon, delta)
    endTimeRU = time()
    
    ourPrecision = len(set(approxFIs) & set(exactFIs)) *1.0 / len(set(approxFIs))
    RUPrecision = len(set(RUapprxFIs) & set(exactFIs)) *1.0 / len(set(RUapprxFIs))
    print('Exact:', endTimeExact-startTimeExact)
    print(ourPrecision, endTimeOurs-startTimeOurs)
    print(RUPrecision, endTimeRU-startTimeRU)
    
'''
Precision test (top-K): select K most frequent pairs. 
Input: dataset file name, # items, top-K, sample increase, delta, alpha
Output: precision, running time.
'''
def testTopKPrecision(fname, I, K, sampleInc, delta, alpha):
    if I*I<K:
        print('Too few items for K')
        return
    singletons = [[i] for i in range(1, I+1)]
    doubletons = []
    if K>=I:
        topKSingle = singletons
    else:
        topKSingle, topKSingleFreqs = exactTopK(fname, singletons, K)
    topKSingle.sort()
    for i in range(min(K,I)):
        for j in range(i+1,min(K,I)):
            doubletons.append(topKSingle[i] + topKSingle[j])
            
    from time import time
    startTimeExact = time()
    exactFIs, exactFreqs = exactTopK(fname, doubletons, K)
    startTimeOurs = endTimeExact = time()
    approxFIs, approxFreqs, nSamples = topK(fname, doubletons, K, sampleInc, delta, alpha)
    startTimeNew = endTimeOurs = time()
    newFIs, newFreqs = newBoundTopK(fname, doubletons, K, .02, .001)
    endTimeNew = time()
    
    ourPrecision = len(set(approxFIs) & set(exactFIs)) *1.0 / len(set(approxFIs))
    newPrecision = len(set(newFIs) & set(exactFIs)) *1.0 / len(set(newFIs))
    print('Exact:', endTimeExact-startTimeExact)
    print(ourPrecision, endTimeOurs-startTimeOurs)
    print(newPrecision, endTimeNew-startTimeNew)
