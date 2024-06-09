# -*- coding: utf-8 -*-
"""
Created on Fri May 17 12:53:28 2024

@author: mustafa.cingiz
"""
import networkx as nx
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
import math
import statistics
import seaborn as sns
import re
#import validationPredAucEnsemble
import matplotlib.pyplot as plt
from sklearn import metrics

diseases=[]
finalEnsembleList=[] # ensemble results, 5 disease, 3 threshold, 15 list
finalSixty=[]   # 60 gene lists, 5 disease, 4 algorithms, 3 threshold


def geneCardEvaluation():
    
    with open('GENECARD-gastric_colorectal_breast_prostate_lung.txt') as f:  #OMIM-EXPANDED-gastric_colorectal_breast_prostate_lung.txt
        for line in f.readlines():
            line=line.strip()
            #print(line)
            splittedLine=line.split(" ")
            diseases.append(splittedLine)
    f.close()  # graph eklendi
    return diseases  #### gene card genes
geneCardEvaluation()
'''
    results=[]
    with open('GENECARD-gastric_colorectal_breast_prostate_lung.txt') as f:  #OMIM-EXPANDED-gastric_colorectal_breast_prostate_lung.txt
        for line in f.readlines():
            line=line.strip()
            print(line)
            splittedLine=line.split(" ")
            diseases.append(splittedLine)
    f.close()  # graph eklendi
    return diseases
'''
################# get results 60 results
sixy_results=pd.read_excel('peerj_all.xlsx', index_col=0)

arr=sixy_results.values.tolist()
resultList= [[k for k in x if str(k) != 'nan'] for x in arr]
################ REGEX
def robustString(myString):
    pattern = r'[^\w\s]'
    cleaned_string = re.sub(pattern, '', myString).strip()
    return cleaned_string.upper()


################# get results 60 results
ensembleResults = pd.read_excel('peerj_ensemble.xlsx', index_col=0)

arrEnsemble = ensembleResults.values.tolist()
##ensemlbeResultList= arrEnsemble[0][3].split(",") #ensemlbeResultList[0][3][1] gibi alınacak değerler

for j in range(15):
    temp = arrEnsemble[j][3].split(",")
    temp=list(map(robustString,temp))
    finalEnsembleList.append(temp)

###########################################  above ensemble preprocess###

for i in range(60):
    ##print(i)
    rltemp=resultList[i][6:]
    temp60=list(map(robustString,rltemp))
    finalSixty.append(temp60)
####################################    all 60 results



prec60=[]
recall60=[]
fmeasure60=[]
tp60=[]
fp60=[]
tn60=[]
fn60=[]
for i in range(60):
    ##print(i)
    rltemp=resultList[i][6:]
    inferredGeneList=list(map(robustString,rltemp))
    if(int(i/12)==0):
        tps=len(list(set(inferredGeneList)& set(diseases[0])))
        fps=len(list(set(inferredGeneList)- set(diseases[0])))
        fns=len(list(set(diseases[0])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        tns=20531- tps-fps-fns
        tpr=recal
        fpr=fps/(tns+fps)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(int(i/12)==1):
        tps=len(list(set(inferredGeneList)& set(diseases[1])))
        fps=len(list(set(inferredGeneList)- set(diseases[1])))
        fns=len(list(set(diseases[1])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        tns=20531- tps-fps-fns
        tpr=recal
        fpr=fps/(tns+fps)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(int(i/12)==2):
        tps=len(list(set(inferredGeneList)& set(diseases[2])))
        fps=len(list(set(inferredGeneList)- set(diseases[2])))
        fns=len(list(set(diseases[2])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        tns=20531- tps-fps-fns
        tpr=recal
        fpr=fps/(tns+fps)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(int(i/12)==3):
        tps=len(list(set(inferredGeneList)& set(diseases[3])))
        fps=len(list(set(inferredGeneList)- set(diseases[3])))
        fns=len(list(set(diseases[3])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        tns=20531- tps-fps-fns
        tpr=recal
        fpr=fps/(tns+fps)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(int(i/12)==4):
        tps=len(list(set(inferredGeneList)& set(diseases[4])))
        fps=len(list(set(inferredGeneList)- set(diseases[4])))
        fns=len(list(set(diseases[4])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        tns=20531- tps-fps-fns
        tpr=recal
        fpr=fps/(tns+fps)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    prec60.insert(i,prec)
    recall60.insert(i,recal)
    fmeasure60.insert(i,fmeasure)
    tp60.insert(i,tps)
    fp60.insert(i,fps)
    tn60.insert(i,tns)
    fn60.insert(i,fns)
    #print("Precision: ", str(prec), " Recall: ", str(recal), "  F-measure:", str(fmeasure), "  TP", str(tps),"  FP:", str(fps),"  FN:", str(fns),"  TN:", str(tns))
dictName = {'precisiom': prec60, 'recall': recall60, 'f-measure': fmeasure60,'tps': tp60, 'fps': fp60, 'fns': fn60,'tns': tn60} 
df60 = pd.DataFrame(dictName)
################################  ROC CURVE VALUES ###########
tpr60=np.array(tp60) / (np.array(tp60)+ np.array(fn60))
fpr60=np.array(tn60) / (np.array(tn60)+ np.array(fp60))

fpr60N=1- fpr60


####################################
precEns=[]
recallEns=[]
fmeasureEns=[]
tpEns=[]
fpEns=[]
tnEns=[]
fnEns=[]
'''
plt.plot(fpr60N[0:3],tpr60[0:3])
fprD=[0]

fprD.extend(fpr60N[0:3])

fprD.append(1)

tprD=[0]

tprD.extend(tpr60[0:3])

tprD.append(1)
plt.plot(fprD,tprD)
##### auc value
my_roc_auc = metrics.auc(fprD,tprD)
'''

##############################ENSEMBLE RESULTS
for i in range(15):
    inferredGeneList=finalEnsembleList[i]
    if( i%5==0):
        tps=len(list(set(inferredGeneList)& set(diseases[0])))
        fps=len(list(set(inferredGeneList)- set(diseases[0])))
        fns=len(list(set(diseases[0])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        tns=20531- tps-fps-fns
        tpr=recal
        fpr=fps/(tns+fps)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(i%5==1):
        tps=len(list(set(inferredGeneList)& set(diseases[1])))
        fps=len(list(set(inferredGeneList)- set(diseases[1])))
        fns=len(list(set(diseases[1])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        tns=20531- tps-fps-fns
        tpr=recal
        fpr=fps/(tns+fps)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(i%5==2):
        tps=len(list(set(inferredGeneList)& set(diseases[2])))
        fps=len(list(set(inferredGeneList)- set(diseases[2])))
        fns=len(list(set(diseases[2])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        tns=20531- tps-fps-fns
        tpr=recal
        fpr=fps/(tns+fps)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(i%5==3):
        tps=len(list(set(inferredGeneList)& set(diseases[3])))
        fps=len(list(set(inferredGeneList)- set(diseases[3])))
        fns=len(list(set(diseases[3])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        tns=20531- tps-fps-fns
        tpr=recal
        fpr=fps/(tns+fps)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(i%5==4):
        tps=len(list(set(inferredGeneList)& set(diseases[4])))
        fps=len(list(set(inferredGeneList)- set(diseases[4])))
        fns=len(list(set(diseases[4])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        tns=20531- tps-fps-fns
        tpr=recal
        fpr=fps/(tns+fps)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    precEns.insert(i,prec)
    recallEns.insert(i,recal)
    fmeasureEns.insert(i,fmeasure)
    tpEns.insert(i,tps)
    fpEns.insert(i,fps)
    tnEns.insert(i,tns)
    fnEns.insert(i,fns)
    #print("Precision: ", str(prec), " Recall: ", str(recal), "  F-measure:", str(fmeasure), "  TP", str(tps),"  FP:", str(fps),"  FN:", str(fns),"  TN:", str(tns))
dictNameEns = {'precisiom': precEns, 'recall': recallEns, 'f-measure': fmeasureEns,'tps': tpEns, 'fps': fpEns, 'fns': fnEns,'tns': tnEns} 
df15 = pd.DataFrame(dictNameEns)
    #print("Precision: ", str(prec), " Recall: ", str(recal), "  F-measure:", str(fmeasure), "  TP", str(tps),"  FP:", str(fps),"  FN:", str(fns),"  TN:", str(tns))

df15.to_excel("ensembleResultsRev1.xlsx")
df60.to_excel("sixtyPredictionResultsRev1.xlsx")

################################  ROC CURVE VALUES ###########
tprEnsemble=np.array(tpEns) / (np.array(tpEns)+ np.array(fnEns))
fprEnsemble=np.array(fpEns) / (np.array(fpEns)+ np.array(tnEns))
##################################  roc curve ###############

plt.plot([0,1],[1,0], 'k--')
plt.plot(tpr60[0], fpr60N[0], label= "gastric_10")
plt.plot(tpr60[1], fpr60N[1], label= "colorectal_10")
plt.plot(tpr60[2], fpr60N[2], label= "breas_10")
plt.plot(tpr60[3], fpr60N[3], label= "prostate_10")
plt.plot(tpr60[4], fpr60N[5], label= "lung_10")

plt.legend()
plt.xlabel("FPR")
plt.ylabel("TPR")
plt.title('Receiver Operating Characteristic')
plt.show()
############################################### ENSEMBLE RESULTS GO PROFILER ##############

precEnsGO=[]
recallEnsGO=[]
fmeasureEnsGO=[]
tpEnsGO=[]
fpEnsGO=[]
tnEnsGO=[]
fnEnsGO=[]

from gprofiler import GProfiler

gp = GProfiler(return_dataframe=True)


for i in range(15):
    inferredGeneList=finalEnsembleList[i]
    if( i%5==0):
        inferredGeneListGO=gp.profile(organism='hsapiens',query=inferredGeneList,sources=["GO:MF","GO:CC","GO:BP"])
        diseasesGO=gp.profile(organism='hsapiens',query=diseases[0],sources=["GO:MF","GO:CC","GO:BP"])
        tpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())& set(diseasesGO.iloc[:,1].values.tolist())))
        fpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())- set(diseasesGO.iloc[:,1].values.tolist())))
        fnsGO=len(list(set(diseasesGO.iloc[:,1].values.tolist())- set(inferredGeneListGO.iloc[:,1].values.tolist())))
        precGO= tpsGO /(tpsGO+fpsGO)
        recalGO= tpsGO/ (tpsGO+fnsGO)
        tnsGO=45000- tpsGO-fpsGO-fnsGO
        tprGO=recalGO
        fprGO=fpsGO/(tnsGO+fpsGO)
        try:
            fmeasureGO= 2*((precGO*recalGO) /(precGO+recalGO))
        except ZeroDivisionError:
            fmeasureG0=0
    elif(i%5==1):
        inferredGeneListGO=gp.profile(organism='hsapiens',query=inferredGeneList,sources=["GO:MF","GO:CC","GO:BP"])
        diseasesGO=gp.profile(organism='hsapiens',query=diseases[1],sources=["GO:MF","GO:CC","GO:BP"])
        tpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())& set(diseasesGO.iloc[:,1].values.tolist())))
        fpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())- set(diseasesGO.iloc[:,1].values.tolist())))
        fnsGO=len(list(set(diseasesGO.iloc[:,1].values.tolist())- set(inferredGeneListGO.iloc[:,1].values.tolist())))
        precGO= tpsGO /(tpsGO+fpsGO)
        recalGO= tpsGO/ (tpsGO+fnsGO)
        tnsGO=45000- tpsGO-fpsGO-fnsGO
        tprGO=recalGO
        fprGO=fpsGO/(tnsGO+fpsGO)
        try:
            fmeasureGO= 2*((precGO*recalGO) /(precGO+recalGO))
        except ZeroDivisionError:
            fmeasureG0=0
    elif(i%5==2):
        inferredGeneListGO=gp.profile(organism='hsapiens',query=inferredGeneList,sources=["GO:MF","GO:CC","GO:BP"])
        diseasesGO=gp.profile(organism='hsapiens',query=diseases[2],sources=["GO:MF","GO:CC","GO:BP"])
        tpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())& set(diseasesGO.iloc[:,1].values.tolist())))
        fpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())- set(diseasesGO.iloc[:,1].values.tolist())))
        fnsGO=len(list(set(diseasesGO.iloc[:,1].values.tolist())- set(inferredGeneListGO.iloc[:,1].values.tolist())))
        precGO= tpsGO /(tpsGO+fpsGO)
        recalGO= tpsGO/ (tpsGO+fnsGO)
        tnsGO=45000- tpsGO-fpsGO-fnsGO
        tprGO=recalGO
        fprGO=fpsGO/(tnsGO+fpsGO)
        try:
            fmeasureGO= 2*((precGO*recalGO) /(precGO+recalGO))
        except ZeroDivisionError:
            fmeasureG0=0
    elif(i%5==3):
        inferredGeneListGO=gp.profile(organism='hsapiens',query=inferredGeneList,sources=["GO:MF","GO:CC","GO:BP"])
        diseasesGO=gp.profile(organism='hsapiens',query=diseases[3],sources=["GO:MF","GO:CC","GO:BP"])
        tpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())& set(diseasesGO.iloc[:,1].values.tolist())))
        fpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())- set(diseasesGO.iloc[:,1].values.tolist())))
        fnsGO=len(list(set(diseasesGO.iloc[:,1].values.tolist())- set(inferredGeneListGO.iloc[:,1].values.tolist())))
        precGO= tpsGO /(tpsGO+fpsGO)
        recalGO= tpsGO/ (tpsGO+fnsGO)
        tnsGO=45000- tpsGO-fpsGO-fnsGO
        tprGO=recalGO
        fprGO=fpsGO/(tnsGO+fpsGO)
        try:
            fmeasureGO= 2*((precGO*recalGO) /(precGO+recalGO))
        except ZeroDivisionError:
            fmeasureG0=0
    elif(i%5==4):
        inferredGeneListGO=gp.profile(organism='hsapiens',query=inferredGeneList,sources=["GO:MF","GO:CC","GO:BP"])
        diseasesGO=gp.profile(organism='hsapiens',query=diseases[4],sources=["GO:MF","GO:CC","GO:BP"])
        tpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())& set(diseasesGO.iloc[:,1].values.tolist())))
        fpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())- set(diseasesGO.iloc[:,1].values.tolist())))
        fnsGO=len(list(set(diseasesGO.iloc[:,1].values.tolist())- set(inferredGeneListGO.iloc[:,1].values.tolist())))
        precGO= tpsGO /(tpsGO+fpsGO)
        recalGO= tpsGO/ (tpsGO+fnsGO)
        tnsGO=45000- tpsGO-fpsGO-fnsGO
        tprGO=recalGO
        fprGO=fpsGO/(tnsGO+fpsGO)
        try:
            fmeasureGO= 2*((precGO*recalGO) /(precGO+recalGO))
        except ZeroDivisionError:
            fmeasureG0=0
    #print("PrecisionG0: ", str(precGO), " RecallG0: ", str(recalGO), "  F-measureG0:", str(fmeasureGO), "  TPG0", str(tpsGO),"  FPG0:", str(fpsGO),"  FNG0:", str(fnsGO),"  TNG0:", str(tnsGO))

    precEnsGO.insert(i,precGO)
    recallEnsGO.insert(i,recalGO)
    fmeasureEnsGO.insert(i,fmeasureGO)
    tpEnsGO.insert(i,tpsGO)
    fpEnsGO.insert(i,fpsGO)
    tnEnsGO.insert(i,tnsGO)
    fnEnsGO.insert(i,fnsGO)
dictNameEnsGO = {'precisiom': precEnsGO, 'recall': recallEnsGO, 'f-measure': fmeasureEnsGO,'tps': tpEnsGO, 'fps': fpEnsGO, 'fns': fnEnsGO,'tns': tnEnsGO} 
df15GO = pd.DataFrame(dictNameEnsGO)
df15GO.to_excel("ensembleGOResultsRev1.xlsx")

################################################################################ 60 GENE LIST GO-TERM OVERLAP ANALYSIS####
'''
precGO=[]
recalGO=[]
fmeasureGO=[]
tpGO=[]
fpGO=[]
tnGO=[]
fnGO=[]
'''
prec60GO=[]
recall60GO=[]
fmeasure60GO=[]
tp60GO=[]
fp60GO=[]
tn60GO=[]
fn60GO=[]
for i in range(60):
    ##print(i)
    rltemp=resultList[i][6:]
    inferredGeneList=list(map(robustString,rltemp))
    if(int(i/12)==0):
        inferredGeneListGO=gp.profile(organism='hsapiens',query=inferredGeneList,sources=["GO:MF","GO:CC","GO:BP"])
        diseasesGO=gp.profile(organism='hsapiens',query=diseases[0],sources=["GO:MF","GO:CC","GO:BP"])
        tpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())& set(diseasesGO.iloc[:,1].values.tolist())))
        fpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())- set(diseasesGO.iloc[:,1].values.tolist())))
        fnsGO=len(list(set(diseasesGO.iloc[:,1].values.tolist())- set(inferredGeneListGO.iloc[:,1].values.tolist())))
        try:
            precGO= tpsGO /(tpsGO+fpsGO)
        except ZeroDivisionError:
            precGO=0
        try:
            recalGO= tpsGO/ (tpsGO+fnsGO)
        except ZeroDivisionError:
            recalGO=0
        tnsGO=45000- tpsGO-fpsGO-fnsGO
        tprGO=recalGO
        fprGO=fpsGO/(tnsGO+fpsGO)
        try:
            fmeasureGO= 2*((precGO*recalGO) /(precGO+recalGO))
        except ZeroDivisionError:
            fmeasureG0=0
    elif(int(i/12)==1):
        inferredGeneListGO=gp.profile(organism='hsapiens',query=inferredGeneList,sources=["GO:MF","GO:CC","GO:BP"])
        diseasesGO=gp.profile(organism='hsapiens',query=diseases[1],sources=["GO:MF","GO:CC","GO:BP"])
        tpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())& set(diseasesGO.iloc[:,1].values.tolist())))
        fpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())- set(diseasesGO.iloc[:,1].values.tolist())))
        fnsGO=len(list(set(diseasesGO.iloc[:,1].values.tolist())- set(inferredGeneListGO.iloc[:,1].values.tolist())))
        try:
            precGO= tpsGO /(tpsGO+fpsGO)
        except ZeroDivisionError:
            precGO=0
        try:
            recalGO= tpsGO/ (tpsGO+fnsGO)
        except ZeroDivisionError:
            recalGO=0
        tnsGO=45000- tpsGO-fpsGO-fnsGO
        tprGO=recalGO
        fprGO=fpsGO/(tnsGO+fpsGO)
        try:
            fmeasureGO= 2*((precGO*recalGO) /(precGO+recalGO))
        except ZeroDivisionError:
            fmeasureG0=0
    elif(int(i/12)==2):
        inferredGeneListGO=gp.profile(organism='hsapiens',query=inferredGeneList,sources=["GO:MF","GO:CC","GO:BP"])
        diseasesGO=gp.profile(organism='hsapiens',query=diseases[2],sources=["GO:MF","GO:CC","GO:BP"])
        tpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())& set(diseasesGO.iloc[:,1].values.tolist())))
        fpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())- set(diseasesGO.iloc[:,1].values.tolist())))
        fnsGO=len(list(set(diseasesGO.iloc[:,1].values.tolist())- set(inferredGeneListGO.iloc[:,1].values.tolist())))
        try:
            precGO= tpsGO /(tpsGO+fpsGO)
        except ZeroDivisionError:
            precGO=0
        try:
            recalGO= tpsGO/ (tpsGO+fnsGO)
        except ZeroDivisionError:
            recalGO=0
        tnsGO=45000- tpsGO-fpsGO-fnsGO
        tprGO=recalGO
        fprGO=fpsGO/(tnsGO+fpsGO)
        try:
            fmeasureGO= 2*((precGO*recalGO) /(precGO+recalGO))
        except ZeroDivisionError:
            fmeasureG0=0
    elif(int(i/12)==3):
        inferredGeneListGO=gp.profile(organism='hsapiens',query=inferredGeneList,sources=["GO:MF","GO:CC","GO:BP"])
        diseasesGO=gp.profile(organism='hsapiens',query=diseases[3],sources=["GO:MF","GO:CC","GO:BP"])
        tpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())& set(diseasesGO.iloc[:,1].values.tolist())))
        fpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())- set(diseasesGO.iloc[:,1].values.tolist())))
        fnsGO=len(list(set(diseasesGO.iloc[:,1].values.tolist())- set(inferredGeneListGO.iloc[:,1].values.tolist())))
        try:
            precGO= tpsGO /(tpsGO+fpsGO)
        except ZeroDivisionError:
            precGO=0
        try:
            recalGO= tpsGO/ (tpsGO+fnsGO)
        except ZeroDivisionError:
            recalGO=0
        tnsGO=45000- tpsGO-fpsGO-fnsGO
        tprGO=recalGO
        fprGO=fpsGO/(tnsGO+fpsGO)
        try:
            fmeasureGO= 2*((precGO*recalGO) /(precGO+recalGO))
        except ZeroDivisionError:
            fmeasureG0=0
    elif(int(i/12)==4):
        inferredGeneListGO=gp.profile(organism='hsapiens',query=inferredGeneList,sources=["GO:MF","GO:CC","GO:BP"])
        diseasesGO=gp.profile(organism='hsapiens',query=diseases[4],sources=["GO:MF","GO:CC","GO:BP"])
        tpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())& set(diseasesGO.iloc[:,1].values.tolist())))
        fpsGO=len(list(set(inferredGeneListGO.iloc[:,1].values.tolist())- set(diseasesGO.iloc[:,1].values.tolist())))
        fnsGO=len(list(set(diseasesGO.iloc[:,1].values.tolist())- set(inferredGeneListGO.iloc[:,1].values.tolist())))
        try:
            precGO= tpsGO /(tpsGO+fpsGO)
        except ZeroDivisionError:
            precGO=0
        try:
            recalGO= tpsGO/ (tpsGO+fnsGO)
        except ZeroDivisionError:
            recalGO=0
        tnsGO=45000- tpsGO-fpsGO-fnsGO
        tprGO=recalGO
        fprGO=fpsGO/(tnsGO+fpsGO)
        try:
            fmeasureGO= 2*((precGO*recalGO) /(precGO+recalGO))
        except ZeroDivisionError:
            fmeasureG0=0
    prec60GO.insert(i,precGO)
    recall60GO.insert(i,recalGO)
    fmeasure60GO.insert(i,fmeasureGO)
    tp60GO.insert(i,tpsGO)
    fp60GO.insert(i,fpsGO)
    tn60GO.insert(i,tnsGO)
    fn60GO.insert(i,fnsGO)
dictName60GO = {'precisiom': prec60GO, 'recall': recall60GO, 'f-measure': fmeasure60GO,'tps': tp60GO, 'fps': fp60GO, 'fns': fn60GO,'tns': tn60GO} 
df60GO = pd.DataFrame(dictName60GO)
df60GO.to_excel("sixtyGOResultsRev1.xlsx")
