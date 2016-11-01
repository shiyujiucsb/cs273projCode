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
    return list(map(lambda x:x*1.0/nTransactions, freqs)), nTransactions

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
        return calcExactFreqs(fname, itemsets)
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
    We do not use sampleInc, epsilon, delta.
Output: top-K itemsets and their frequencies
'''
def exactTopK(fname, itemsets, K, sampleInc, epsilon, delta):
    freqs, n = calcExactFreqs(fname, itemsets)
    topKFIs = sorted(range(len(itemsets)), key=lambda x:-freqs[x])
    return list(map(lambda x:(itemsets[x]), topKFIs))[:K], sorted(freqs, key=lambda x:-x)[:K], n

'''
Approx top-k algorithm with new bound.
Input: dataset file name, itemsets, top-K, epsilon, delta
    We do not use sampleInc.
Output: top-K itemsets and their frequencies
'''
def newBoundTopK(fname, itemsets, K, sampleInc, epsilon, delta):
    freqs, n = newBoundApprox(fname, itemsets, epsilon, delta)
    topKFIs = sorted(range(len(itemsets)), key=lambda x:-freqs[x])
    return list(map(lambda x:(itemsets[x]), topKFIs))[:K], sorted(freqs, key=lambda x:-x)[:K], n

'''
Approx top-k algorithm by RU's (epsilon, delta)-approximation method.
Input: dataset file name, itemsets, top-K, epsilon, delta
    We do not use sampleInc.
Output: top-K itemsets and their frequencies
'''
def RUTopK(fname, itemsets, K, sampleInc, epsilon, delta):
    freqs, n = RUApprox(fname, itemsets, epsilon, delta)
    topKFIs = sorted(range(len(itemsets)), key=lambda x:-freqs[x])
    return list(map(lambda x:(itemsets[x]), topKFIs))[:K], sorted(freqs, key=lambda x:-x)[:K], n

'''
Our progressive top-K approximate algorithm
Input: dataset file name, itemsets, top-K, sample increase, epsilon, delta
Output: top-K itemsets and their frequencies, # samples
'''
def ProgressiveTopK(fname, itemsets, K, sampleInc, epsilon, delta):
    freqs = [0 for _ in range(len(itemsets))]
    f = open(fname, 'r')
    lines = f.readlines()
    nTransactions = len(lines)
    f.close()
    from random import randint
    from math import log, exp
    N=min(nTransactions,int((log(2*len(itemsets))-log(delta))/(2*epsilon*epsilon))+1)
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
        if errorProb <= delta or n>N:
            return list(map(lambda x:itemsets[x], \
                            sorted(range(len(itemsets)), key=lambda x:-freqs[x])[:K])),\
                    list(map(lambda x:x*1.0/n, tmp[:K])), n

'''
RU's classical top-K approximate algorithm
Input: dataset file name, itemsets, top-K, sample increase, epsilon, delta
Output: top-K itemsets and their frequencies, # samples
'''
def RUProgressiveTopK(fname, itemsets, K, sampleInc, epsilon, delta):
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
        ell = max(freqs)**.5
        tmp = sorted(freqs, key=lambda x:-x)
        middle = (tmp[K-1]+tmp[K])*.5
        errorProb=0.0
        for i in range(len(tmp)):
            dProb = exp(-2.0/n*abs(tmp[i]-middle)**2)
            errorProb +=dProb
            if i>=K:
                if errorProb>delta or errorProb+(len(tmp)-1-i)*dProb<=delta:
                    break
        Delta = 2*ell/n*(2*log(len(itemsets)))**0.5 + (2*log(2/delta)/n)**0.5
        if errorProb <= delta or n>nTransactions or Delta<=epsilon:
            return list(map(lambda x:itemsets[x], \
                            sorted(range(len(itemsets)), key=lambda x:-freqs[x])[:K])),\
                    list(map(lambda x:x*1.0/n, tmp[:K])), n

'''
A-Priori Top-K algorithm
Input: dataset file name, # items, top-K, function to compute Tok-K, sample increase, epsilon, delta
Output: top-K itemsets and their frequencies, # samples
'''
def AprioriTopK(fname, I, K, topKFunction, sampleInc, epsilon, delta):
    if I*I<K:
        print('Too few items for K')
        return
    singletons = [[i] for i in range(1, I+1)]
    doubletons = []
    if K>=I:
        topKSingle = singletons
    else:
        topKSingle, topKSingleFreqs, nSmp = topKFunction(fname, singletons, K, sampleInc, epsilon, delta)
    topKSingle.sort()
    for i in range(min(K,I)):
        for j in range(i+1, min(K,I)):
            doubletons.append(topKSingle[i] + topKSingle[j])
    topKDouble, topKDoubleFreqs, nSmp = topKFunction(fname, doubletons, K, sampleInc, epsilon, delta)
    return topKDouble, topKDoubleFreqs, nSmp

'''
Precision/recall test (non-progressive top-K): select K most frequent pairs. 
Input: dataset file name, # items, top-K, sample increase, epsilon, delta
Output: precision, running time.
'''
def testNonProgressiveTopK(fname, I, K, sampleInc, epsilon, delta):
    from time import time
    startTimeExact = time()
    exactFIs, exactFreqs, n = AprioriTopK(fname, I, K, exactTopK, sampleInc, epsilon, delta)
    startTimeOurs = endTimeExact = time()
    approxFIs, approxFreqs, nUs = AprioriTopK(fname, I, K, newBoundTopK, sampleInc, epsilon, delta)
    startTimeRU = endTimeOurs = time()
    RUapprxFIs, RUapprxFreqs, nRU = AprioriTopK(fname, I, K, RUTopK, sampleInc, epsilon, delta)
    endTimeRU = time()

    setApproxFIs = set([str(i) for i in approxFIs])
    setExactFIs = set([str(i) for i in exactFIs])
    setRUapprxFIs = set([str(i) for i in RUapprxFIs])
    
    ourPrecision = len(setApproxFIs & setExactFIs) *1.0 / len(setApproxFIs)
    RUPrecision = len(setRUapprxFIs & setExactFIs) *1.0 / len(setRUapprxFIs)
    print(fname)
    #print('Exact:', endTimeExact-startTimeExact)
    print(ourPrecision, endTimeOurs-startTimeOurs, nUs)
    print(RUPrecision, endTimeRU-startTimeRU, nRU)
    print('')
    
'''
Precision/recall test (progressive top-K): select K most frequent pairs. 
Input: dataset file name, # items, top-K, sample increase, epsilon, delta
Output: precision, running time.
'''
def testProgressiveTopK(fname, I, K, sampleInc, epsilon, delta):
    from time import time
    startTimeExact = time()
    exactFIs, exactFreqs, n = AprioriTopK(fname, I, K, exactTopK, sampleInc, epsilon, delta)
    startTimeOurs = endTimeExact = time()
    approxFIs, approxFreqs, nUs = AprioriTopK(fname, I, K, ProgressiveTopK, sampleInc, epsilon, delta)
    startTimeRU = endTimeOurs = time()
    RUapprxFIs, RUapprxFreqs, nRU = AprioriTopK(fname, I, K, RUProgressiveTopK, sampleInc, epsilon, delta)
    endTimeRU = time()

    setApproxFIs = set([str(i) for i in approxFIs])
    setExactFIs = set([str(i) for i in exactFIs])
    setRUapprxFIs = set([str(i) for i in RUapprxFIs])
    
    ourPrecision = len(setApproxFIs & setExactFIs) *1.0 / len(setApproxFIs)
    RUPrecision = len(setRUapprxFIs & setExactFIs) *1.0 / len(setRUapprxFIs)
    print(fname)
    #print('Exact:', endTimeExact-startTimeExact)
    print(ourPrecision, endTimeOurs-startTimeOurs, nUs)
    print(RUPrecision, endTimeRU-startTimeRU, nRU)

# test zone

#testNonProgressiveTopK('dataset/accidents.dat.txt',68,100,100,0.05,0.0001)
#testNonProgressiveTopK('dataset/chess.dat.txt',75,100,100,0.05,0.0001)
testNonProgressiveTopK('dataset/connect.dat.txt',70,100,100,0.05,0.0001)
#testNonProgressiveTopK('dataset/kosarak.dat.txt',70,100,100,0.05,0.0001)
#testNonProgressiveTopK('dataset/mushroom.dat.txt',70,100,100,0.05,0.0001)
#testNonProgressiveTopK('dataset/pumsb.dat.txt',70,100,100,0.05,0.0001)
#testNonProgressiveTopK('dataset/pumsb_star.dat.txt',70,100,100,0.05,0.0001)
#testNonProgressiveTopK('dataset/retail.dat.txt',70,100,100,0.05,0.0001)

print("========")

#testProgressiveTopK('dataset/accidents.dat.txt',68,100,100,0.05,0.0001)
#testProgressiveTopK('dataset/chess.dat.txt',75,100,100,0.05,0.0001)
#testProgressiveTopK('dataset/connect.dat.txt',70,100,100,0.05,0.0001)
#testProgressiveTopK('dataset/kosarak.dat.txt',70,100,100,0.05,0.0001)
#testProgressiveTopK('dataset/mushroom.dat.txt',70,100,100,0.05,0.0001)
#testProgressiveTopK('dataset/pumsb.dat.txt',70,100,100,0.05,0.0001)
#testProgressiveTopK('dataset/pumsb_star.dat.txt',70,100,100,0.05,0.0001)
#testProgressiveTopK('dataset/retail.dat.txt',70,100,100,0.05,0.0001)
