# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 20:36:14 2023

@author: özgür-pc
"""
import networkx as nx
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
import math
import statistics
import seaborn as sns


'''
file = open('hprd-csv.csv')
hprdData = csv.reader(file)
G = nx.Graph()

rows = []
for row in hprdData:
    rows.append(row)
    G.add_edge(row[1], row[2])
rows

#hprdNumpy = np.asarray(rows)

#G.add_nodes_from(rows[:,1], rows[:,2])

'''
G = nx.Graph()
genler=[]
with open('hprd_unique_final.txt') as f:
    for line in f.readlines():
        splittedLine=line.split("\t")
        firstGene=splittedLine[1].strip()
        genler.append(firstGene)
        secondGene=splittedLine[2].strip()
        genler.append(secondGene)
        G.add_edge(firstGene, secondGene)
f.close()  # graph eklendi

ugenler=np.array(genler)
uniqueGenes=np.unique(ugenler) # 12890 tekrarsiz genler bulundu
#print(len(uniqueGenes))
G.add_nodes_from(list(uniqueGenes)) # dugumler graph'a eklendi

def linkPred(method, disease, topK):
    
    if(method=="pac"):
        preds = nx.preferential_attachment(G)
    elif(method=="adamicAdar"):
        preds = nx.adamic_adar_index(G)
    elif(method=="jaccard"):
        preds = nx.jaccard_coefficient(G)
    elif(method=="common_neighbour"):
        preds = nx.resource_allocation_index(G)
       
    #dct = {'firstGene':'Peter', 'secondGene':'male', 'pa_score':21}  # gastric_cancer_base_list
    dct={}
    list_pac=[]
    
    if(disease=="gastric"):
        base_gene_list=["IL1RN","IL1B,""IRF1","KLF6","APC","PIK3CA","CASP10","CDH1","ERBB2","MUTYH","FGFR2"]
    elif(disease=="colorectal"):
        base_gene_list= ["CRCS6","CRCS7","DCC","MLH1","CRCS5","CRCS2","PLA2G2A","AXIN2","BRAF","MSH2","SMAD7","MLH3","MCC","TGFBR2","CRCS8","PDGFRL","CRCS9","TLR2","HMPS1","APC","PIK3CA","DLC1","BUB1B","MSH6","BAX","TLR4","FGFR3","TP53","CTNNB1","CRCS11","CRCS10","FLCN","NRAS","CCND1","EP300","CHEK2","AKT1","BUB1","PMS2","PMS1"]
    elif(disease=="breast"):
        base_gene_list= ["NQO1","XRCC3","BRCD2","PALB2","ESR1","RAD51A","PPM1D","PIK3CA","ATM","SLC22A1L","TP53","KRAS2","TSG101","PHB","HMMR","RB1CC1","BRCA2","BRCA3","BRIP1","CASP8","RAD54L","CDH1","CHEK2","BRCD1","AKT1","BRCATA","BCPR","BARD1"]
    elif(disease=="prostate"):
        base_gene_list=["MSR1","ZFHX3","ELAC2","HPCQTL19","AR","EHBP1","KLF6","HPC3","HPC5","HPC4","HPC7","HPC6","HPC9","MAD1L1","HIP1","HPC11","HPC10","RNASEL","PCAP","CD82","HPC15","HPC14","HPCX2","PTEN","BRCA2","HPCX1","MXI1","CHEK2","EPHB2","MSMB"]
    elif(disease=="lung"):
        base_gene_list=["TSG11","CHRNA3","CHRNA5","DDX26","MPO","BRAF","DLEC1","LNCR1","EGFR","RASSF1","IRF1","LNCR4","CASP8","LNCR3","CYP2A6","PIK3CA","PPP2R1B","SLC22A1L","ERCC6","MAP3K8","KRAS2"]
    
    
    for u, v, p in preds:
        if(u in base_gene_list or v in base_gene_list):
            dct={}
            dct["firstGene"]=u
            dct["secondGene"]=v
            dct["pa_score"]=p
            list_pac.append(dct)
            #print(f"({u}, {v}) -> {p}")
    
    pac_scores = [list_pac[i]["pa_score"] for i in range(len(list_pac))]
    #sns.boxplot(pac_scores,orient="v") # SEABORN, orient v veya h olacak
    #sns.stripplot(pd.DataFrame(pac_scores),marker="o", alpha=0.3, color="black")
    pac_scores.sort()
    
    #pac_score_median=statistics.median(new_list)
    #pac_score_threshold= pac_scores[round(len(pac_scores)*0.8)] # RATIO BAZLI bu deger 0.8 degil parametrik yapacagiz
    topK=-1*topK
    pac_score_threshold= pac_scores[topK] # top 10'daki degeri aliyor
    #pac_score_threshold= pac_scores[-10] # top 10'daki degeri aliyor
    
    
    # BU HASTALIKLA ILISKILI TUUMMMM GENLERIN HEPSININ AA, JC vb. DEGERINE BAKARAK HEPSININ ICINDEN EN YUKSEK 10 SKORU OLANI ALIYOR
    imp_int_PAC=[list_pac[i] for i in range(len(list_pac)) if list_pac[i]["pa_score"]>=pac_score_threshold]
    
    topGenes=[]
    for i in range(len(imp_int_PAC)):
        topGenes.append(imp_int_PAC[i]["firstGene"])
        topGenes.append(imp_int_PAC[i]["secondGene"])
        
    topDiseaseRelatedGenes=list(set(topGenes)- set(base_gene_list))
    topDiseaseRelatedGenes= list(set(topDiseaseRelatedGenes))
    return topDiseaseRelatedGenes


def validationPred(inferredGeneList, disease):
    diseases=[]
    with open('gastric_colorectal_breast_prostate_lung.txt') as f:
        for line in f.readlines():
            line=line.strip()
            splittedLine=line.split("\t")
            diseases.append(splittedLine)
    f.close()  # graph eklendi
        
    if(disease=="gastric"):
        tps=len(list(set(inferredGeneList)& set(diseases[0])))
        fps=len(list(set(inferredGeneList)- set(diseases[0])))
        fns=len(list(set(diseases[0])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(disease=="colorectal"):
        tps=len(list(set(inferredGeneList)& set(diseases[1])))
        fps=len(list(set(inferredGeneList)- set(diseases[1])))
        fns=len(list(set(diseases[1])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(disease=="breast"):
        tps=len(list(set(inferredGeneList)& set(diseases[2])))
        fps=len(list(set(inferredGeneList)- set(diseases[2])))
        fns=len(list(set(diseases[2])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(disease=="prostate"):
        tps=len(list(set(inferredGeneList)& set(diseases[3])))
        fps=len(list(set(inferredGeneList)- set(diseases[3])))
        fns=len(list(set(diseases[3])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    elif(disease=="lung"):
        tps=len(list(set(inferredGeneList)& set(diseases[4])))
        fps=len(list(set(inferredGeneList)- set(diseases[4])))
        fns=len(list(set(diseases[4])- set(inferredGeneList)))
        prec= tps /(tps+fps)
        recal= tps/ (tps+fns)
        try:
            fmeasure= 2*((prec*recal) /(prec+recal))
        except ZeroDivisionError:
            fmeasure=0
    print("Precision: ", str(prec), " Recall: ", str(recal), "  F-measure:", str(fmeasure))
    return prec,recal,fmeasure, inferredGeneList





#predictions= linkPred(method="adamicAdar", disease="colorectal", topK=20)
#validationPred(predictions, disease="colorectal")


algos=["pac","adamicAdar","jaccard","common_neighbour"]
disorders=["gastric","colorectal","breast","prostate","lung"]
allKTops=[10,50,100]

allResults=[]
final_parameters=[]

for i in disorders:
    for j in algos:
        for k in allKTops:
            predictions= linkPred(method=j, disease=i, topK=k)
            result=validationPred(predictions, disease=i)
            allResults.append(result)
            final_parameters.append([i,j,k])

#x=[1,2.3,4,["aaa","bbb","ccc"],5,6.2,11,["acd","ert"]]
counter=1
with open(r'parametersOto.txt', 'w') as fp:
    for item in allResults:
        
        # write each item on a new line
        
        #if(counter % 5==0):
        #    fp.write("\n")
        fp.write(str(item)+ " ")
        fp.write("\n")
        counter=counter+1
    fp.close()            

import pickle
with open('final_parameters.pkl', 'wb') as frr:
    pickle.dump(final_parameters, frr)

###################################################
import pandas as pd

mydata=pd.read_excel('label_data.xlsx') 

#### top10

topten= mydata[::3]
topfifty=mydata[1::3]
tophundred=mydata[2::3]

# clean nan values
cleanedList = [x for x in topten.iloc[0,7:].tolist() if str(x) != 'nan']
# remove special characters [,],"" and  space
#cleanString= cleanedList[0].replace("'", "").replace("[","").replace("]","").replace(" ","")
def eleminateSpecChar(mystr):
    mystr= mystr.replace("'", "").replace("[","").replace("]","").replace(" ","")
    return mystr

output_str = list(map(eleminateSpecChar, cleanedList))

### TOPTEN ICIN ####

allgenesTopTen=[]

counter=0
for i in range(5):
    tempL=[]
    for j in range(4):
        cleanedList = [x for x in topten.iloc[counter,7:].tolist() if str(x) != 'nan']
        output_str = list(map(eleminateSpecChar, cleanedList))
        for k in range(len(output_str)):
            tempL.append(output_str[k])
        counter=counter+1
    allgenesTopTen.append(tempL)

# dictionary 

gastric_dct10={x:list(allgenesTopTen[0]).count(x) for x in list(allgenesTopTen[0]) }
colorectal_dct10={x:list(allgenesTopTen[1]).count(x) for x in list(allgenesTopTen[1]) }
breast_dct10={x:list(allgenesTopTen[2]).count(x) for x in list(allgenesTopTen[2]) }
prostate_dct10={x:list(allgenesTopTen[3]).count(x) for x in list(allgenesTopTen[3]) }
lung_dct10={x:list(allgenesTopTen[4]).count(x) for x in list(allgenesTopTen[4]) }






allgenesTopFifty=[]

counter=0
for i in range(5):
    tempL=[]
    for j in range(4):
        cleanedList = [x for x in topfifty.iloc[counter,7:].tolist() if str(x) != 'nan']
        output_str = list(map(eleminateSpecChar, cleanedList))
        for k in range(len(output_str)):
            tempL.append(output_str[k])
        counter=counter+1
    allgenesTopFifty.append(tempL)

# dictionary 
tem=set(allgenesTopFifty[0])
gastric_dct50={x:list(allgenesTopFifty[0]).count(x) for x in list(allgenesTopFifty[0]) }
colorectal_dct50={x:list(allgenesTopFifty[1]).count(x) for x in list(allgenesTopFifty[1]) }
breast_dct50={x:list(allgenesTopFifty[2]).count(x) for x in list(allgenesTopFifty[2]) }
prostate_dct50={x:list(allgenesTopFifty[3]).count(x) for x in list(allgenesTopFifty[3]) }
lung_dct50={x:list(allgenesTopFifty[4]).count(x) for x in list(allgenesTopFifty[4]) }




allgenesTopHundred=[]

counter=0
for i in range(5):
    tempL=[]
    for j in range(4):
        cleanedList = [x for x in tophundred.iloc[counter,7:].tolist() if str(x) != 'nan']
        output_str = list(map(eleminateSpecChar, cleanedList))
        for k in range(len(output_str)):
            tempL.append(output_str[k])
        counter=counter+1
    allgenesTopHundred.append(tempL)

# dictionary 

gastric_dct100={x:list(allgenesTopHundred[0]).count(x) for x in list(allgenesTopHundred[0]) }
colorectal_dct100={x:list(allgenesTopHundred[1]).count(x) for x in list(allgenesTopHundred[1]) }
breast_dct100={x:list(allgenesTopHundred[2]).count(x) for x in list(allgenesTopHundred[2]) }
prostate_dct100={x:list(allgenesTopHundred[3]).count(x) for x in list(allgenesTopHundred[3]) }
lung_dct100={x:list(allgenesTopHundred[4]).count(x) for x in list(allgenesTopHundred[4]) }

wholeRanked=[]
wholeRanked.append(gastric_dct10)
wholeRanked.append(colorectal_dct10)
wholeRanked.append(breast_dct10)
wholeRanked.append(prostate_dct10)
wholeRanked.append(lung_dct10)
wholeRanked.append(gastric_dct50)
wholeRanked.append(colorectal_dct50)
wholeRanked.append(breast_dct50)
wholeRanked.append(prostate_dct50)
wholeRanked.append(lung_dct50)
wholeRanked.append(gastric_dct100)
wholeRanked.append(colorectal_dct100)
wholeRanked.append(breast_dct100)
wholeRanked.append(prostate_dct100)
wholeRanked.append(lung_dct100)

###########################
ensm_Results=[]   ## ALL ENSEMBLE RESULTS
myc=0

for ind in range(len(wholeRanked)):
    myc=myc+1
    tempS=[]
    for u,v in wholeRanked[ind].items():#for u,v in gastric_dct100.items():
        if v>1:
            tempS.append(u)
              
    tempS=list(set(tempS))
        
    if(myc%5==1):
        secondPar="gastric"
    elif(myc%5==2):
        secondPar="colorectal"
    elif(myc%5==3):
        secondPar="breast"
    elif(myc%5==4):
        secondPar="prostate"
    elif(myc%5==0):
        secondPar="lung"
    prec,recal,fmeasure, inferredGeneList= validationPred(tempS,secondPar)
    addInfo=[ind,secondPar,prec,recal,fmeasure,inferredGeneList]
    ensm_Results.append(addInfo)
    


## ALL ENSEMBLE RESULTS WRITE EXCEL
dfensm_Results = pd.DataFrame(ensm_Results)
writer = pd.ExcelWriter('ensembleResults.xlsx', engine='xlsxwriter')
dfensm_Results.to_excel(writer, sheet_name='welcome', index=False)
writer.save()   
    
    

#################top write top10 to excel
dfgastric_dct10 = pd.DataFrame(data=gastric_dct10, index=["Gene- Count"])
dfgastric_dct10=dfgastric_dct10.T
dfgastric_dct10.to_excel("dfgastric_dct10.xlsx", index=True)

dfcolorectal_dct10 = pd.DataFrame(data=colorectal_dct10, index=["Gene- Count"])
dfcolorectal_dct10=dfcolorectal_dct10.T
dfcolorectal_dct10.to_excel("dfcolorectal_dct10.xlsx", index=True)

dfbreast_dct10 = pd.DataFrame(data=breast_dct10,index=["Gene- Count"])
dfbreast_dct10=dfbreast_dct10.T
dfbreast_dct10.to_excel("dfbreast_dct10.xlsx", index=True)

dfprostate_dct10 = pd.DataFrame(data=prostate_dct10,index=["Gene- Count"])
dfprostate_dct10=dfprostate_dct10.T
dfprostate_dct10.to_excel("dfprostate_dct10.xlsx", index=True)

dflung_dct10 = pd.DataFrame(data=lung_dct10,index=["Gene- Count"])
dflung_dct10=dflung_dct10.T
dflung_dct10.to_excel("dflung_dct10.xlsx", index=True)


#################top write top50 to excel
dfgastric_dct50 = pd.DataFrame(data=gastric_dct50, index=["Gene- Count"])
dfgastric_dct50=dfgastric_dct50.T
dfgastric_dct50.to_excel("dfgastric_dct50.xlsx", index=True)

dfcolorectal_dct50 = pd.DataFrame(data=colorectal_dct50, index=["Gene- Count"])
dfcolorectal_dct50=dfcolorectal_dct50.T
dfcolorectal_dct50.to_excel("dfcolorectal_dct50.xlsx", index=True)

dfbreast_dct50 = pd.DataFrame(data=breast_dct50,index=["Gene- Count"])
dfbreast_dct50=dfbreast_dct50.T
dfbreast_dct50.to_excel("dfbreast_dct50.xlsx", index=True)

dfprostate_dct50 = pd.DataFrame(data=prostate_dct50,index=["Gene- Count"])
dfprostate_dct50=dfprostate_dct50.T
dfprostate_dct50.to_excel("dfprostate_dct50.xlsx", index=True)

dflung_dct50 = pd.DataFrame(data=lung_dct50,index=["Gene- Count"])
dflung_dct50=dflung_dct50.T
dflung_dct50.to_excel("dflung_dct50.xlsx", index=True)


#################top write top100 to excel
dfgastric_dct100 = pd.DataFrame(data=gastric_dct100, index=["Gene- Count"])
dfgastric_dct100=dfgastric_dct100.T
dfgastric_dct100.to_excel("dfgastric_dct100.xlsx", index=True)

dfcolorectal_dct100 = pd.DataFrame(data=colorectal_dct100, index=["Gene- Count"])
dfcolorectal_dct100=dfcolorectal_dct100.T
dfcolorectal_dct100.to_excel("dfcolorectal_dct100.xlsx", index=True)

dfbreast_dct100 = pd.DataFrame(data=breast_dct100,index=["Gene- Count"])
dfbreast_dct100=dfbreast_dct100.T
dfbreast_dct100.to_excel("dfbreast_dct100.xlsx", index=True)

dfprostate_dct100 = pd.DataFrame(data=prostate_dct100,index=["Gene- Count"])
dfprostate_dct100=dfprostate_dct100.T
dfprostate_dct100.to_excel("dfprostate_dct100.xlsx", index=True)

dflung_dct100 = pd.DataFrame(data=lung_dct100,index=["Gene- Count"])
dflung_dct100=dflung_dct100.T
dflung_dct100.to_excel("dflung_dct100.xlsx", index=True)

#convert into excel


#["gastric","colorectal","breast","prostate","lung"]
#def cleanFromNan(snan):
    


#output_str = list(map(eleminateSpecChar, cleanedList))



#cleanedList = [x for x in topten.iloc[0,7:].tolist() if str(x) != 'nan']


## number of row and column list
#rows = len(group1)
#columns = len(group1[0])








        
