#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:42:06 2020

main script for ODT optimisation
needs DataProcessing.py to function

@author: martijnkorf
"""

from docplex.mp.model import Model
import DataProcessing
import numpy

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------

dataset = "student"         # choose from: mush, ttt, heloc, votes, student, bc, heart, monks-1, krkp
treeTopology = "depth3"     # choose from: depth2 depth2.5 depth3 imbalanced
trainingVSClass = 300       # choose from: 90%/600 or an integer value
specialOption = ""          # choose from: nothing, allTrees, tryTechniques, subsampling, crossValidation

# tree model options
relaxation = True
strengthening = True
anchoring = True
numericalFeature = True
imbalanceCor = False
maxRank = False             # alternative sampling technique
corrHeatmap = False         # print feature selection information

# fetch and configure model parameters and tree parameters
features, featuresNamed, anchorFeatures, ysamples, J, G, N, aMatrix, yourResults, GNumerical, lastjperg, datasetSparsity, C, Iplus, Imin, Ntrain = DataProcessing.readData(dataset, trainingVSClass, maxRank, imbalanceCor)
B, K, Bplus, Bmin, Kleft, Kright, KnoLeaves, Kleaves, Bright, Bleft, treeTopology = DataProcessing.Topology(treeTopology)
print("Configurations: %s %s %s" % (dataset, treeTopology, str(trainingVSClass)))

# ----------------------------------------------------------------------------
# Main Functions for the DT Model
# ----------------------------------------------------------------------------

def build_Tree(relaxation, strengthening, anchoring, numericalFeature):
    mdl = Model("Decision Tree %s" % dataset)
    print("constructing model with %r %r %r %r..." % 
          (relaxation,strengthening,anchoring, numericalFeature))
    
    # create the variables
    z = [[0 for k in range(K)] for J in range(J)]
    if (relaxation) :
        z1 = mdl.continuous_var_matrix(J,Kleaves,name="z")
        z2 = mdl.binary_var_matrix(J,KnoLeaves,name="z")
        for j in range(J) :
            for k in Kleaves   : z[j][k] = z1[(j,k)]
            for k in KnoLeaves : z[j][k] = z2[(j,k)]
        v = mdl.continuous_var_matrix(G, K,0,1, name="v")
        c = mdl.continuous_var_matrix(N, B,0,1, name="c")
    else :
        z3 = mdl.binary_var_matrix(J, K, name="z")
        for j in range(J) :
            for k in range(K) : z[j][k] = z3[j,k]
        v = mdl.binary_var_matrix(G, K, name="v")
        c = mdl.binary_var_matrix(N, B, name="c")
    L,R = [[0 for k in range(K)] for i in range(Ntrain)],[[0 for k in range(K)] for i in range(Ntrain)]
    for k in range(K) : 
        for i in range(Ntrain) :
            L[i][k] = sum(aMatrix[i][j]*z[j][k] for j in range(J))
            R[i][k] = 1 - L[i][k]
    print("created all variables", end = '')
    
    # add the restrictions of S
    mdl.add_constraints(mdl.sum(v[g,k] for g in range(G)) == 1 for k in range(K))
    mdl.add_constraints(z[j][k] <= v[(featuresNamed[j].group-1),k] for k in range(K) for j in range(J))
    
    # add restrictions of Q
    if (strengthening) :
        mdl.add_constraints(mdl.sum(c[(i,b)] for b in Bleft[k])  <= L[i][k] for i in range(Ntrain) for k in range(K))
        mdl.add_constraints(mdl.sum(c[(i,b)] for b in Bright[k]) <= R[i][k] for i in range(Ntrain) for k in range(K))
    else :
        mdl.add_constraints(c[(i,b)] <= L[i][k] for b in range(B) for i in range(Ntrain) for k in Kleft[b])
        mdl.add_constraints(c[(i,b)] <= R[i][k] for b in range(B) for i in range(Ntrain) for k in Kright[b])
    mdl.add_constraints(mdl.sum(c[(i,b)] for b in range(B)) == 1 for i in range(Ntrain))
    print(", all restrictions")
    
    #add anchoring
    if (anchoring == True) and (treeTopology == "depth2" or treeTopology == "depth3") :
        mdl.add_constraints(z[anchorFeatures[g]][k] == v[(g,k)] for g in range(G) for k in KnoLeaves)
    elif (anchoring == True) :
        if   (treeTopology == "depth2.5")   : Ksym = {1}
        elif (treeTopology == "imbalanced") : Ksym = {2}
        mdl.add_constraints(z[anchorFeatures[g]][k] == v[(g,k)] for g in range(G) for k in Ksym)
    
    # add numerical constraint, if any
    if (len(GNumerical) > 0 and numericalFeature) :
        listNum = [i for i in GNumerical]
        print("numerical restraints added for %s ; C = %f" % (listNum, C))
        w = mdl.binary_var_matrix(K, GNumerical, name="w")
        mdl.add_constraints(z[j][k] >= (z[jplus1][k] - w[(k,g)]) for k in range(K) for g in GNumerical 
                            for j in range(lastjperg[g]-9,lastjperg[g])
                            for jplus1 in range(lastjperg[g]-8,lastjperg[g]+1))
        mdl.add_constraints(z[j][k] >= (z[jmin1][k] - (1 - w[(k,g)])) for k in range(K) for g in GNumerical
                            for j in range(lastjperg[g]-8,lastjperg[g]+1)
                            for jmin1  in range(lastjperg[g]-9,lastjperg[g]))
    else : print("only categorical data ; C = %f" % C)
    print("Topo: %s ; Bmin: %s ; Bplus: %s" % (treeTopology, Bmin, Bplus))
    # Maximize objective function, here maybe flag to which one to optimize :)    
    objective_function = mdl.sum( c[i,b] for i in Iplus for b in Bplus) + C * mdl.sum( c[i,b] for i in Imin for b in Bmin)
    mdl.maximize(objective_function)

    return mdl

# prints a solved tree
def print_Tree() :
    
    print("The Decision Tree (%s):" %(dataset))
    #Print (node,group,features)
    for node in range(K) :
        features = [] ; groupBranch = "no group found" #reset to find errors
        for group in range(G) : 
            vkgString = "v_%i_%i" % (group,node)
            if (curSol[vkgString] == 1 or (curSol[vkgString] - 1 < 0.01 and curSol[vkgString]>0)):
                groupBranch = group
                break
        features = [feature for feature in range(J) 
                    if curSol["z_%i_%i" % (feature,node)] == 1] # store selected feature
        print("node: %i ; group: %2i ; features: %s"
              % (node,groupBranch,str(features)[1:-1]))
        
    training_accuracy = (curMdl.objective_value / N ) * 100
    print("Objective Value: %4.2d %% ; Solve Time: %4.5f ; Nodes: %i"
          % (training_accuracy,curSol.solve_details.time, curSol.solve_details.nb_nodes_processed))
  
# calculates the classification acc of train and test, returns test accuracy          
def classification() :
    
    sizeSamples = [N, len(aMatrix) - N]
    boundsSample = [0, sizeSamples[0], sizeSamples[0] + sizeSamples[1]]
    text = ["Train. Accuracy:","Class. Accuracy:"]
    treeTopologies = ["depth2","depth2.5","depth3","imbalanced"]
    treeTopo = treeTopologies.index(treeTopology)
    nextNode = [[[2,"0","0"],[1,"1","1"]],
                [[4,3,"0","0","0"],[1,2,"1","1","1"]],
                [[4,3,"0","0",6,"0","0"],[1,2,"1","1",5,"1","1"]],
                [[6,5,4,"0","0","0","0"],[1,2,3,"1","1","1","1"]]]
    
    #first check train. acc then calculate class. acc
    for time in range(2) :   
        counterCorrect,counterError = 0,0
        for i in range(boundsSample[time], boundsSample[time+1]):
            currentNode = 0 #reset node
            for turn in range(4) : # 4 is max depth of tree
                
                left = sum(aMatrix[i][feature] * curSol["z_%i_%i" % (feature, currentNode)]
                           for feature in range(J))
                # rounding errors
                if   ((1 - left) < 0.01) : left = 1
                elif ((0 + left) < 0.01) : left = 0
                if (left != 0) and (left != 1): print("invalid left: %f node: %i"
                                                    % (left, currentNode))
                currentNode = nextNode[treeTopo][left][currentNode]
                if (currentNode == "1" or currentNode == "0") :
                    break
                
            if (ysamples[i] == int(currentNode)): counterCorrect += 1
            else : counterError += 1
        
        accuracy = ( counterCorrect / sizeSamples[time] ) * 100
        print("%s %4.2f %% ; size: %4.i ; # false %3.i "
              % (text[time], accuracy, sizeSamples[time], counterError))
        #if (time == 0) : trainAcc = accuracy
    print("")
    
    return accuracy

# ----------------------------------------------------------------------------
# Calling all the scripts
# ----------------------------------------------------------------------------

# tNumerical Feature and correlation Imbalance tables code  
if (specialOption == "allTreesB"):
    print("Solving for all tree sizes")
    
    topologies = ["depth2","depth2.5","depth3","imbalanced"]
    numFeat = [True,False]
    numTop  = len(topologies)
    corrImbalance = [True,False]
    two = len(numFeat)
    repeats = 2
    solveTs,nodesProc = numpy.zeros((repeats,two,len(topologies))), numpy.zeros((repeats,two,len(topologies)))
    objValues,classAc = numpy.zeros((repeats,two,len(topologies))), numpy.zeros((repeats,two,len(topologies)))
    
    for rep in range(repeats) :
        sampleSize1 = "90%/600" # 90% of 600 pls fix
        print("sample:", sampleSize1)
        features, featuresNamed, anchorFeatures, ysamples, J, G, aMatrix, yourResults, GNumerical, lastjperg, datasetSparsity, C, Iplus, Imin, Ntrain = DataProcessing.readData(dataset, trainingVSClass, maxRank, imbalanceCor)
        N = Ntrain ; 
        if (imbalanceCor == False): C = 1
        
        for num in range(two) :

            for t in range(numTop) :
                
                B, K, Bplus, Bmin, Kleft, Kright, KnoLeaves, Kleaves, Bright, Bleft, treeTopology = DataProcessing.Topology(treeTopology[t])
                curMdl = build_Tree(True,True,True,True)
                curMdl.print_information()
                print( "solving model for %s with num %r ..." % (topologies[t],numFeat[num]))
                curMdl.set_time_limit(600) 
                curSol = curMdl.solve()
                nodesProc[rep][num][t] =curSol.solve_details.nb_nodes_processed
                #objValues[rep][num][t] = (curMdl.objective_value / Ntrain ) * 100
                print_Tree()
                classAc[rep][num][t],objValues[rep][num][t] = classification()
                curMdl.clear()

    avgTimes   = [[(sum(solveTs[i][num][tree]    for i in range(repeats)) / repeats) 
                   for tree in range(numTop)] for num in range(two)]
    avgNodes   = [[(sum(nodesProc[i][num][tree]  for i in range(repeats)) / repeats) 
                   for tree in range(numTop)] for num in range(two)]
    avgObj     = [[(sum(objValues[i][num][tree]  for i in range(repeats)) / repeats) 
                   for tree in range(numTop)] for num in range(two)]
    avgClassAc = [[(sum(classAc[i][num][tree]    for i in range(repeats)) / repeats) 
                   for tree in range(numTop)] for num in range(two)]
    print(avgTimes,avgNodes)
    print(avgObj,avgClassAc)

# code that will estimate models with different improvement techniques
elif (specialOption == "tryTechniques") :
    
    print("Solving for different improvements")
    techniques = [[False,False,False],[True,True,False],[False,True,True],[True,False,True],[True,True,True]]
    numTech = len(techniques)
    repeats = 3
    solveTs,nodesProc,objValues,classAc = numpy.zeros((repeats,numTech)), numpy.zeros((repeats,numTech)),numpy.zeros((repeats,numTech)),numpy.zeros((repeats,numTech))
    
    for rep in range(repeats) :
        features, featuresNamed, anchorFeatures, ysamples, J, G, aMatrix, yourResults, GNumerical, lastjperg, datasetSparsity, C, Iplus, Imin, Ntrain = DataProcessing.readData(dataset, trainingVSClass, maxRank, imbalanceCor)

        for tech in range(numTech) :
        
            print("Solving for: %r %r %r" % (techniques[tech][0], 
                                             techniques[tech][1], techniques[tech][2]))
            curMdl = build_Tree(techniques[tech][0],techniques[tech][1],
                                techniques[tech][2],False)
            curMdl.print_information()
            print( "solving model ...\n")
            curMdl.set_time_limit(600)
            curSol = curMdl.solve()
            solveTs[rep][tech] = curSol.solve_details.time
            nodesProc[rep][tech] = curSol.solve_details.nb_nodes_processed
            objValues[rep][tech] = (curMdl.objective_value / Ntrain ) * 100
            print_Tree()
            classAc[rep][tech]   = classification()
            
    avgTimes   = [sum(solveTs[i][tech]   for i in range(repeats)) / repeats
                  for tech in range(numTech)]
    avgNodes   = [sum(nodesProc[i][tech] for i in range(repeats)) / repeats 
                  for tech in range(numTech)]
    avgObj     = [sum(objValues[i][tech] for i in range(repeats)) / repeats
                  for tech in range(numTech)]
    avgClassAc = [sum(classAc[i][tech]   for i in range(repeats)) / repeats
                  for tech in range(numTech)]
    print(avgTimes,avgNodes,avgObj,avgClassAc)
    
# code that will perform K-Fold crossValidation, adjust numK for different # of folds
elif (specialOption == "crossValidation") :
    from sklearn.model_selection import StratifiedKFold
    
    features, featuresNamed, anchorFeatures, ysamplesOG, J, G, N, aMatrixOG, yourResults, GNumerical, lastjperg, datasetSparsity, C, Iplus, Imin, Ntrain = DataProcessing.readData(dataset, trainingVSClass, maxRank,imbalanceCor)
    numK = 5 ; k = 0
    kfold = StratifiedKFold(numK,False)
    solveTs,nodesProc,objValues,classAc = [],[],[],[]
    print("Solving with K-fold Crossvalidation, K=%i" % numK)
    
    for train,test in kfold.split(aMatrixOG,ysamplesOG):
        print("train size: %i ; test size: %i" % (len(train), len(test)))
        #create new sorted Amatrix,ysamples,Iplus,Imin
        aMatrix  = [aMatrixOG[i] for i in train]
        ysamples = [ysamplesOG[i] for i in train]
        for i in test:
            aMatrix.append(aMatrixOG[i])
            ysamples.append(ysamplesOG[i])
        Iplus = [i for i in range(len(train)) if ysamples[i] == 1]
        Imin =  [i for i in range(len(train)) if ysamples[i] == 0]
        Ntrain = len(train)
        # build, solve and print tree
        curMdl = build_Tree(True,True,True,True)
        curMdl.print_information()
        print( "solving model ...\n")
        curMdl.set_time_limit(600)
        curSol = curMdl.solve()
        solveTs.append(curSol.solve_details.time)
        nodesProc.append(curSol.solve_details.nb_nodes_processed)
        objValues.append( (curMdl.objective_value / Ntrain ) * 100 )
        print_Tree()
        classAc.append( classification() )
    
    avgTimes   = [sum(solveTs[i]   for i in range(numK)) / numK]
    avgNodes   = [sum(nodesProc[i] for i in range(numK)) / numK]
    avgObj     = [sum(objValues[i] for i in range(numK)) / numK]
    avgClassAc = [sum(classAc[i]   for i in range(numK)) / numK]
    print(avgTimes,avgNodes,avgObj,avgClassAc)
    
# code that will perform random subsampling, to compare with cross validation
elif (specialOption == "randomSubsampling") :
    
    repeats = 5
    Ntrain = int( 0.8*len(yourResults) )
    Ntest = len(yourResults) - Ntrain
    solveTs,nodesProc,objValues,classAc = [],[],[],[]
    print("Solving with random subsampling, repeats=%i" % repeats)
    
    for rep in range(repeats):
        # fetch new random sample
        features, featuresNamed, anchorFeatures, ysamples, J, G, N, aMatrix, yourResults, GNumerical, lastjperg, datasetSparsity, _, Iplus, Imin, Ntrain = DataProcessing.readData(dataset, Ntrain, maxRank, imbalanceCor)
        print("train size: %i ; test size: %i" % (Ntrain, Ntest))
        # build, solve and print tree
        curMdl = build_Tree(True,True,True,True)
        curMdl.print_information()
        print( "solving model ...\n")
        curMdl.set_time_limit(600)
        curSol = curMdl.solve()
        solveTs.append(curSol.solve_details.time)
        nodesProc.append(curSol.solve_details.nb_nodes_processed)
        objValues.append( (curMdl.objective_value / Ntrain ) * 100 )
        print_Tree()
        classAc.append( classification() )
    
    avgTimes   = [sum(solveTs[i]   for i in range(repeats)) / repeats]
    avgNodes   = [sum(nodesProc[i] for i in range(repeats)) / repeats]
    avgObj     = [sum(objValues[i] for i in range(repeats)) / repeats]
    avgClassAc = [sum(classAc[i]   for i in range(repeats)) / repeats]
    print(avgTimes,avgNodes,avgObj,avgClassAc)
    
# estimate a single model with certain settings
else :
    curMdl = build_Tree(relaxation,strengthening,anchoring,numericalFeature)
    curMdl.print_information()
    print( "solving model ...\n")
    curMdl.set_time_limit(600)
    curSol = curMdl.solve()
    print_Tree()
    classification()

# create diagrams for feature selection
if (corrHeatmap) :
    
    import seaborn as sns
    import pandas  as pd
    from sklearn.feature_selection import SelectKBest
    from sklearn.feature_selection import chi2
    from sklearn.ensemble import ExtraTreesClassifier
    import matplotlib.pyplot as plt
    
    data1 = aMatrix.copy()
    for i in range(len(aMatrix)):
        data1[i].append(ysamples[i])
    data = pd.DataFrame(data1)
    X = data.iloc[:,0:J]  #independent columns
    y = data.iloc[:,-1]    #target column i.e price range
    
    model = ExtraTreesClassifier()
    model.fit(X,y)
    #plot graph of feature importances for better visualization
    feat_importances = pd.Series(model.feature_importances_, index=X.columns)
    feat_importances.nlargest(min(25,J)).plot(kind='barh')
    plt.show()
    
    #get correlations of each features in dataset
    corrmat = data.corr()
    top_corr_features = corrmat.index
    plt.figure(figsize=(20,20))
    #plot heat map
    g=sns.heatmap(data[top_corr_features].corr(),cmap="RdYlGn",vmin=-1,vmax=1,annot=True)
    
    bestfeatures = SelectKBest(score_func=chi2, k=min(25,J))
    fit = bestfeatures.fit(X,y)
    dfscores = pd.DataFrame(fit.scores_)
    dfcolumns = pd.DataFrame(X.columns)
    #concat two dataframes for better visualization 
    featureScores = pd.concat([dfcolumns,dfscores],axis=1)
    featureScores.columns = ['Specs','Score']  #naming  the dataframe columns
    #print(featureScores.nlargest(J,'Score'))  #print 10 best features
    