# -*- coding: utf-8 -*-

import math
import random
import numpy
import matplotlib.pyplot as plt
import copy
import csv
import pickle
import sys
import matplotlib.ticker as tkr
#import seaborn as sns
import pandas as pd
#sns.set(context="paper", font="monospace")

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

def solutions():
    sol = []
    gene1 = [1,1,0,1]
    gene2 = [0,0,0,0]
    gene3 = [0,0,0,1]
    gene4 = [0,0,1,1]
    gene5 = [0,0,1,0]

    sol.append(gene1)
    sol.append(gene2)
    sol.append(gene3)
    sol.append(gene4)
    sol.append(gene5)

    gene1fit = 1.0
    gene2fit = 0.5
    gene3fit = 0.5
    gene4fit = 0.6
    gene5fit = 0.4

    return sol

def gene(geneSize):
    gene = []
    for i in range (0, geneSize):
        gene.append(random.randint(0,1));
    return gene

def initPop(popSize, geneSize):
    pop = []
    for i in range (0,popSize): 
        pop.append(gene(geneSize))
    return pop

def fixedGene(geneSize):
    gene = []
    for i in range (0, geneSize):
        if (i < (geneSize * 0.5)-0.01):
            gene.append(0)
        else: 
            gene.append(1)
    return gene

def initFixedPop(popSize, geneSize):
    pop = []
    for i in range (0,popSize): 
        pop.append(fixedGene(geneSize))
    return pop

def mutateGene(gene, mutRate):
    for i in range(0, len(gene)):
        r = random.uniform(0,1)
        if (r < mutRate):
            gene[i] = random.randint(0,1);
    return gene

def mutatePop(pop, popfit, mutRate):
    mutPop = []
    mutPopFit = []
    for i in range(0, len(pop)):
        randomInd = random.randint(0, len(pop)-1)
        if popfit[randomInd] > -0.1:
            gene = copy.deepcopy(pop[randomInd])
            mutGene = mutateGene(gene, mutRate)
            if (crossover == True):
                if (random.randrange(0.0,1.0) < crossoverRate):
                    otherInd = random.randint(0,len(pop)-1)
                    if (popfit[otherInd] > -0.1 and otherInd != randomInd):
                        mutGene = twoPointCrossover(gene, pop[otherInd])
            mutPop.append(mutGene)
            mutPopFit.append(evaluateGene(mutGene))
        else:
            mutGene = copy.deepcopy(pop[randomInd])
            mutPop.append(copy.deepcopy(pop[randomInd]))
            mutPopFit.append(copy.deepcopy(popfit[randomInd]))
    return mutPop, mutPopFit
    
def evaluateGene(gene):
    # H-IFF
    fitness = 0;
    # First the amount of layers in the hierarchy is determined. 
    layers = 0
    for i in range(10):
        i += 1
        if pow(2, i) == len(gene):
            layers = i
            break
    if layers == 0:
        print("ERROR, layer must be 2 to power something for H-IFF algorithm to work")
    currentLayer = 0

    # calculating fitness (Note that in the original H-IFF algorithm the 
    while currentLayer < layers:
        currentLayer += 1
        compareArray = []
        bitlength = pow(2, (currentLayer -1))
        i = 0
        while i < (len(gene)-bitlength + 1):
            cBitString = []
            for j in range(bitlength):
                cBitString.append(gene[j+i])
            i += bitlength
            compareArray.append(cBitString)
        amountChecks = len(compareArray) / 2
        i = 0
        while i < (len(compareArray)-1):
            i += 1
            if (compareArray[i] == compareArray[i-1]):
                fitness += pow(2,(currentLayer-1))
            i += 1
    return fitness

def getdistr(popfit, mf):
    colorRuns = plt.cm.rainbow(numpy.linspace(0, 1, 20))
    distr = [0] * len(colorRuns)
    for i in range(len(popfit)):
        length = len(colorRuns) -1
        type = int(round(float(popfit[i]) / float(mf)  * float(length)))
        distr[type] += 1
    return distr

def evaluatePop(pop):
    popFit = []
    for i in range(0, len(pop)):
        popFit.append(evaluateGene(pop[i]))
    return popFit

def comparePop(pop, popfit, mutatedPop, mutatedpopfit):
    for i in range(0,len(pop)):
        r = random.randint(0,len(mutatedPop)-1)
        if (popfit[r] < mutatedpopfit[i]):
           popfit[r] = mutatedpopfit[i]
           pop[r] = mutatedPop[i]
    return pop, popfit

def agepop(popage, popfit, maxage):
    for i in range(0, len(popage)):
        popage[i] += 1
 #       if (popage[i] > maxage):
 #           popfit[i] = -1
    return popage, popfit

def agepopDP(popage, popfit, dp):
    for i in range(0, len(popage)):
        r = (random.uniform(0,1))
        if (r < dp):
            popfit[i] = -1.0
    return popfit

def initpopage(pop):
    popage = []
    for i in range(0, len(pop)):
        popage.append(0)
    return popage

def run(ps, gs, mr, ag, dp, ar, mf):
    populationSize = ps
    geneSize = gs
    mutRate = mr
    amountGenerations = ag
    maxage = 100
    deathprob = dp
    amountruns = ar 
    allaverages = []
    allMax = []
    distribution = []
    avav = []
    avMax = []
    minEr = []
    maxEr = []
    maxMinEr = []
    maxMaxEr = []
    pop = []
    popFit =[]
    
    allPops = []
    allPopFits = []
    for n in range (0,amountruns):
        curCompletion = 0;
        print "run", n
        pop = initPop(populationSize, geneSize)
        popfit = []
        popfit = evaluatePop(pop)
        pops = []
        popfits = []
        popage = []
        popage = initpopage(pop)
        av = numpy.average(popfit)
        alldistr = []
        alldistr.append(getdistr(popfit, mf))
        averages = []
        maxFits = []

        maxFits.append(numpy.max(av))
        averages.append(av)
        for i in range (0,amountGenerations): 
            if (i % extinctionInterval == 0 and i != 0):
                popfit = extinctionEvent(popage,popfit)
      #      print(str((100.0/float(amountGenerations)*float(i))) + " complete")
            if ((float(100) / float(amountGenerations) * float(i)) > curCompletion):
                curCompletion += 1
        #        print(str(curCompletion) + "% complete")
            popfit = agepopDP(popage,popfit,deathprob)
            mutatedPop = []
            mutPopFit = []
            mutatedPop, mutPopFit = mutatePop(pop, popfit, mutRate)
            pop, popfit = comparePop(pop, popfit, mutatedPop, mutPopFit)

            av = numpy.average(popfit)
            averages.append(av)

            maxFit = 0
            for j in range(len(popfit)):
                if (popfit[j] > maxFit):
                    maxFit = popfit[j]
            maxFits.append(maxFit)
            distr = []
            if (i % interval == 0):
                pops.append(copy.deepcopy(pop))
                popfits.append(copy.deepcopy(popfit))
                distr = getdistr(popfit, mf)
                alldistr.append(distr)
        allaverages.append(averages)    
        allMax.append(maxFits)
        distribution.append(alldistr)
        averages = []
        pops.append(copy.deepcopy(pop))
        popfits.append(copy.deepcopy(popfit))  
        allPops.append(copy.deepcopy(pops))
        allPopFits.append(copy.deepcopy(popfits))
       
    
    for n in range(0,amountGenerations + 1):
        sum = 0
        erArray = []
        for i in range(amountruns):
            sum += allaverages[i][n]
            erArray.append(allaverages[i][n])
        avr = sum/amountruns
        avav.append(avr)
        minEr.append(numpy.percentile(erArray,25))
        maxEr.append(numpy.percentile(erArray,75))
        if ((25 / amountGenerations * i) + 50 > curCompletion):
            curCompletion += 1
            print(curCompletion, "% complete")

    for n in range(0,amountGenerations + 1):
        sum = 0    
        maxErArray = []
        for i in range(amountruns):
            sum += allMax[i][n]
            maxErArray.append(allMax[i][n])
        avr = sum/amountruns
        avMax.append(avr)
        maxMinEr.append(numpy.percentile(maxErArray,25))
        maxMaxEr.append(numpy.percentile(maxErArray,75))
        if (((25 / amountGenerations * i) + 75) > curCompletion):
            curCompletion += 1
            print(curCompletion, "% complete")
    
    return allaverages, avav, allMax, avMax, minEr, maxEr, maxMinEr, maxMaxEr, distribution, allPops, allPopFits

def plotDistribution():
    fig2 = plt.figure()
    ax1 = fig2.add_subplot(111)

def randrange(n, vmin, vmax): #temp
    '''
    Helper function to make an array of random numbers having shape (n, )
    with each number distributed Uniform(vmin, vmax).
    '''
    return (vmax - vmin)*numpy.random.rand(n) + vmin

def plotGraph(allaverages, avav, allMax, avMax, minEr, maxEr, maxMinEr, maxMaxEr, distr, dp, interval, maxFit, finalPops, finalPopFits, titleText):    
    colorRuns = plt.cm.rainbow(numpy.linspace(0, 1, len(allMax)))

    fig = plt.figure()
    ax1 = fig.add_subplot(221)   #top left
    for i in range(len(allaverages)):
        plt.plot(allaverages[i], color = colorRuns[i])
    plt.plot(avav, color='black', linewidth=3)
    plt.plot(minEr, linestyle=':', color = 'black', linewidth=3)
    plt.plot(maxEr, linestyle=':', color = 'black', linewidth=3)
    
    plt.ylabel("fitness")
    plt.xlabel("generations")
    ax1.title.set_text(titleText)
    axes = plt.gca()
    axes.set_ylim([0,maxFit+1])

    ax11 = fig.add_subplot(223)
    for i in range(len(allMax)):
        plt.plot(allMax[i], color = colorRuns[i])
    plt.plot(avMax, color='black', linewidth=3)
    plt.plot(maxMinEr, linestyle=':', color = 'black', linewidth=3)
    plt.plot(maxMaxEr, linestyle=':', color = 'black', linewidth=3)
    
    plt.ylabel("fitness")
    plt.xlabel("generations")
    plt.title("max")

    ax1.title.set_text(titleText)
    axes = plt.gca()
    axes.set_ylim([0,maxFit + 1])
    
    indTypesColor = plt.cm.rainbow(numpy.linspace(0, 1, len(distr[0][0])))
    types = []
    for i in range(0, len(distr[0])):
        tmp_av = []
        for n in range(len(indTypesColor)):
            sum = 0
            for j in range(len(distr)):
                sum += distr[j][i][n]
            tmp_av.append(sum)
        types.append(tmp_av)

    ax2 = fig.add_subplot(222)

    plt.xlim(0,int(len(types)*interval))
    plt.title("Population Fitness Distribution")
    for n in range(len(types)):
        lengthOfSomething = len(types[n])
        currentBottom = 0
 #       ind = numpy.arange(lengthOfSomething)
        ar = numpy.array(types[n])  
        for i in range(len(types[n])):
            plt.bar(int(n * interval),ar[i], 0.8 * float(interval),color = indTypesColor[i], bottom = currentBottom, edgecolor = None)
            currentBottom+=numpy.array(types[n][i])
    ax3 = fig.add_subplot(224)
    plt.xlim(-len(finalPops[0][0][0]) -1, len(finalPops[0][0][0]) + 1)
    plt.ylim(0,maxFit)
    for p in range(len(finalPops)):
        popNum = len(finalPops[0])
        distr_1 = []
        distr_2 = []
        for n in range(len(finalPops[p])):
            distr_1 = []
            for i in range(len(finalPops[p][n])):
                amountOnes = 0
                amountZeros = 0
                #print finalPops[n][i]
                for j in range(len(finalPops[p][n][i])):
                    if(finalPops[p][n][i][j] == 1):
                        amountOnes += 1
                    else:
                        amountZeros += 1
#        if (amountOnes > amountZeros):
                distr_1.append(amountOnes-amountZeros)
           # distr_2.append(copy.deepcopy(distr_1))
            alp = 0.0    
            if (n != 0):
               alp = float(n) / float(popNum) 
            plt.scatter(distr_1, finalPopFits[p][n], s = 100.0, alpha=alp, c = colorRuns[p])

def maxFoundFun(allMax, distr, maxFit, finalPops, finalPopFits, titleText, filename):
    interval = len(allMax[0]) / 40;
    print len(allMax[0]) 
    print 'plotting simple graph'
    colorRuns = plt.cm.rainbow(numpy.linspace(0, 1, len(allMax)))
    foundMaxFit = []
    maxFits = []
    for i in range(20):
        foundMaxFit.append(False)
        maxFits.append(0.0)
    amountMaxFit = 0
    
    max100Er = []
    genArray = []
    amountRuns = len(allMax)
    gen = len(allMax[0])
    for i in range(gen):
        sum = 0
        a = []
        for j in range(amountRuns):
            sum += allMax[j][i]
            av = float(sum/amountRuns)
            a.append(allMax[j][i])
            if (allMax[j][i] == maxFit):
                foundMaxFit[j] = True
            if (allMax[j][i] > maxFits[j]):
                maxFits[j] = allMax[j][i]
        max100Er.append(numpy.percentile(a,100))
        genArray.append(i)

    for i in range(20):
        if (foundMaxFit[i] == True):
            amountMaxFit += 1
        #print("Run ", i, " : ", foundMaxFit[i])
    print "maxfits: " + str(maxFits)
    print "amount found solution: " , amountMaxFit
  


def plotGraphSimple(allMax, distr, maxFit, finalPops, finalPopFits, titleText):
    interval = len(allMax[0]) / 40;
    print len(allMax[0]) 
    print 'plotting simple graph'
    colorRuns = plt.cm.rainbow(numpy.linspace(0, 1, len(allMax)))
    foundMaxFit = []
    for i in range(200):
        foundMaxFit.append(False)

    fig = plt.figure()
    ax11 = fig.add_subplot(111)
    yArr = [0,16,32,48,64,80,96,112,128]
    minor_ticks =  numpy.arange(0,(maxFit+16),1)
    minor_ticks_x = numpy.arange(0,(len(allMax[0])+1),int((len(allMax[0])*0.02)))
    major_ticks = numpy.arange(0,(maxFit+16),16)
    ax11.set_yticks(major_ticks)
    ax11.set_yticks(minor_ticks, minor = True)
    ax11.set_xticks(minor_ticks_x, minor = True)
    avMax = []
    #min25Er = []
    #min0Er = []
    #max75Er = []
    #max100Er = []
    amountRuns = len(allMax)
    gen = len(allMax[0])
    for i in range(gen):
        sum = 0
        for j in range(amountRuns):
            sum += allMax[j][i]
            av = float(sum/amountRuns)
            if (allMax[j][i] == maxFit):
                foundMaxFit[j] = True;
        avMax.append(av)
        #min25Er.append(numpy.percentile(allMax[j]),25)
        #min0Er.append(numpy.percentile(allMax[j]),0)
        #max75Er.append(numpy.percentile(allMax[j]),75)
        #max100Er.append(numpy.percentile(allMax[j]),100)
    for i in range(len(allMax)):
        plt.plot(allMax[i], color = colorRuns[i])
    ax11.set_ylim(0, maxFit + 1)
    plt.rc('grid', linestyle="-", color='black')
    #plt.plot(mir0Er)
    #plt.plot(max100Er)

    #ax11.set_yticks(major_ticks)
    #ax11.tick_params(axis = 'both', which = 'minor', labelsize = 0)

    #ax11.grid(b=True, which='major', alpha = 0.6, linestyle='-')
    #ax11.grid(which='minor', alpha = 0.3, linestyle='-')
    
    plt.plot(avMax, color='black', linewidth=3)

#    plt.plot(maxMinEr, linestyle=':', color = 'black', linewidth=3)
#    plt.plot(maxMaxEr, linestyle=':', color = 'black', linewidth=3)
    for tick in ax11.xaxis.get_major_ticks():
                tick.label.set_fontsize(16) 
    for tick in ax11.yaxis.get_major_ticks():
                tick.label.set_fontsize(16) 
    plt.ylabel("fitness", fontsize = 16)
    plt.xlabel("generation", fontsize = 16)
#    plt.title("max")

#    ax11.title.set_text(titleText)
    axes = plt.gca()
#    axes.set_ylim([0,maxFit + 1])

    for i in range(200):
        print("Run ", i, " : ", foundMaxFit[i])

    print 'plotted max '
    
    indTypesColor = plt.cm.rainbow(numpy.linspace(0, 1, len(distr[0][0])))
    types = []
    for i in range(0, len(distr[0])):
        tmp_av = []
        for n in range(len(indTypesColor)):
            sum = 0
            for j in range(len(distr)):
                sum += distr[j][i][n]
            tmp_av.append(sum)
        types.append(tmp_av)
    
    fig2 = plt.figure();
    ax2 = fig2.add_subplot(211)

    plt.xlim(0,int(len(types)*interval))
    #plt.title("Population Fitness Distribution")
    for n in range(len(types)):
        lengthOfSomething = len(types[n])
        currentBottom = 0
        ar = numpy.array(types[n])  
        for i in range(len(types[n])):
            plt.bar((n * interval),ar[i], 0.8 * float(interval),color = indTypesColor[i], bottom = currentBottom, edgecolor = None)
            currentBottom+=numpy.array(types[n][i])

    print 'plotted distr'

    ax3 = fig2.add_subplot(212)
    plt.xlim(-1, (len(finalPops[0][0][0])) + 1)
    ax3.xaxis.set_ticks(numpy.arange(0, (len(finalPops[0][0][0])) +1, 8))
    plt.ylim(0,maxFit + 1)
    major_ticks = numpy.arange(0,(maxFit+16),32)
    ax3.set_yticks(major_ticks)
    
    for p in range(len(finalPops)):
        popNum = len(finalPops[0])
        distr_1 = []
        distr_2 = []
        for n in range(len(finalPops[p])):
            distr_1 = []
            for i in range(len(finalPops[p][n])):
                amountOnes = 0
                amountZeros = 0
                for j in range(len(finalPops[p][n][i])):
                    if(finalPops[p][n][i][j] == 1):
                        amountOnes += 1
                    else:
                        amountZeros += 1
                distr_1.append(amountOnes)
            alp = 0.0    
            if (n != 0):
               alp = float(n) / float(popNum) 
            plt.scatter(distr_1, finalPopFits[p][n], s = 100.0, alpha=alp, c = colorRuns[p])
    print 'plotted distr'

    fig3D = plt.figure();
    ax3d = fig3D.add_subplot(111,projection='3d')
    ax3d.set_ylim(-1, (len(finalPops[0][0][0]))+1)
    ax3d.yaxis.set_ticks(numpy.arange(0, (len(finalPops[0][0][0])) + 1, 16))
    for tick in ax3.xaxis.get_major_ticks():
        tick.label.set_fontsize(18) 
    for tick in ax3.yaxis.get_major_ticks():
        tick.label.set_fontsize(18) 
    ax3.set_xlabel('distribution', fontsize = 18)
    ax3.set_ylabel('fitness', fontsize = 18)


    ax3d.set_xlim(0,len(finalPops[0])*interval / 1000);
    ax3d.set_zlim(0,maxFit + 1)
    for p in range(len(finalPopFits)):
        for n in range(len(finalPops[p])):
            length = []
            fitnessChumps = []
            distr_1 = []
            for i in range(len(finalPops[p][n])):
                amountOnes = 0
                amountZeros = 0
                for j in range(len(finalPops[p][n][i])):
                    if(finalPops[p][n][i][j] == 1):
                        amountOnes += 1
                    else:
                        amountZeros += 1
                distr_1.append(amountOnes)
                length.append(n * interval / 1000)
            ax3d.scatter(length, distr_1, finalPopFits[p][n] , zdir='z',c = colorRuns[p], marker='o')
    
    major_ticks = numpy.arange(0,(maxFit+16),32)
    ax3d.set_zticks(major_ticks)
    ax3d.set_xlabel('generation * 10^3',fontsize = 14)
    ax3d.set_ylabel('distribution',fontsize = 14)
    ax3d.set_zlabel('fitness',fontsize = 14)
    for tick in ax3d.xaxis.get_major_ticks():
        tick.label.set_fontsize(14) 
    for tick in ax3d.yaxis.get_major_ticks():
        tick.label.set_fontsize(14) 
    for tick in ax3d.zaxis.get_major_ticks():
        tick.label.set_fontsize(14) 
   
def main():
    populationSize = 50
    geneSize = 32
    mutRate = 0.1
    amountGenerations = 5000
    deathprob = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]
    amountruns = 5
    global interval
    if amountGenerations < 20:
        interval = 1
    else:
        interval = int(amountGenerations*0.025)
    global extinctionInterval
    global extinctionRate
    global crossover
    global crossoverRate
    crossover = False
    crossoverRate = 0.000
    extinctionInterval = 1000000
    extinctionRate = 0.0
    extinctionIntervalArray = [10] #cannot be 0
    global optimizationType 
	
    avav = []
    allaverages = []
    minEr = []
    maxEr = []

    print sys.argv
    print len(sys.argv)
    print 'loading arguments if any'    

    fittestGeneHiff = gene(geneSize)
    for i in range(geneSize):
        fittestGeneHiff.append(0)        
    maxFit = evaluateGene(fittestGeneHiff)
            
    if (1 == 0):
        for n in range(65):
            switch = 0;
            fittestGeneHiff = [];
            for i in range(geneSize):
                fittestGeneHiff.append(0);
            for f in range(n):
                fittestGeneHiff[(int(f * 2) + f) % 64] = 1
            maxFit = evaluateGene(fittestGeneHiff)
            print maxFit, n#, fittestGeneHiff
            amountZeros = 0
            for i in range(geneSize):
                if fittestGeneHiff[i] == 0:
                    amountZeros += 1
    for m in range(len(deathprob)):
        allaverages, avav, allMax, avMax, minEr, maxEr, maxMinEr, maxMaxEr, distr, finalPops, finalPopFits = run(populationSize, geneSize, mutRate, amountGenerations, deathprob[m], amountruns, maxFit)
        name = 'DR ' + str(deathprob[m])
        plotGraph(allaverages, avav, allMax, avMax, minEr, maxEr, maxMinEr, maxMaxEr, distr, 0.0, interval, maxFit, finalPops, finalPopFits, name)
        #plotGraphSimple(allMax, distr, maxFit, finalPops, finalPopFits, name)
    plt.show()

    sol = solutions()

main()
