# -*- coding: utf-8 -*-
"""
Created on Sat May 18 10:05:26 2024

@author: mustafa.cingiz
"""

import numpy as np
import pandas as pd
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_curve, auc

# go-term evaluation
mydata = pd.read_excel('sixtyGOResultsComplete.xlsx')

# genecard overlap analysis
#mydata = pd.read_excel('sixtyPredictionResultsComplete.xlsx')

#mydata = pd.read_excel('sixtyPredictionResultsComplete.xlsx') # 
myarray=mydata.iloc[:,:].values
#temp= myarray[0,4]  # ornek veri 1 tane

tps=myarray[:,6]
fps=myarray[:,7]
fns=myarray[:,8]
tns=myarray[:,9]
tpr=tps/ (tps+fns)
fpr=fps/(tns+fps)

tprL=[]
fprL=[]
diseases=["GC","CRC","BRC","PRC","LC"]
algorithms=["PAC","AA","JACCARD","RAI","SMV"]

cntr=0
for i in range(0,75,3):
    tprL.append(0)
    tprL.extend(tpr[i:i+3])
    tprL.append(1)
    fprL.append(0)
    fprL.extend(fpr[i:i+3])
    fprL.append(1)

tprL=list(map(float, tprL))
fprL=list(map(float, fprL))

# Alt grafikler için figür ve eksenleri oluştur
fig, axs = plt.subplots(3, 2, figsize=(14, 18)) #fig, axs = plt.subplots(3, 2, figsize=(14, 18)) ---_> 3 satir 2 sutunken 

# Her alt grafik için ROC eğrilerini çiz
for i in range(5):
    ax = axs[i // 2, i % 2]  #  ax = axs[i // 3, i % 3] #ax = axs[i // 2, i % 2] ---_> 3 satir 2 sutunken  
    
    #for #name, model in models.items():
    algoInd=0
    for j in range(5): #for j in range(0,20,5):
        
        #model.fit(X_train, y_train)
        #y_score = model.predict_proba(X_test)[:, 1]
        #fpr, tpr, _ = roc_curve(y_test, y_score)
        #roc_auc = auc(fpr, tpr)
        
        #ax.plot(fpr, tpr, lw=2, label=f'{name} (area = {roc_auc:.2f})')  items.sort(reverse=True)
        try:
            #roc_auc = auc(fprL[i*5+j:i*5+j+5], tprL[i*5+j:i*5+j+5])
            roc_auc = auc(fprL[i*25+j*5:i*25+j*5+5], tprL[i*25+j*5:i*25+j*5+5])
            roc_score="{:.3f}".format(roc_auc)
    
        except:
            continue
        try:
            ax.plot(fprL[i*25+j*5:i*25+j*5+5], tprL[i*25+j*5:i*25+j*5+5], lw=2, label= '{},{}, auc_score= {}'.format(diseases[i],algorithms[j],roc_score))
        except:
            continue
        algoInd=algoInd+1
      
    ax.plot([0, 1], [0, 1], color='gray', lw=3, linestyle='--')
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title(' ROC Curve of {}'.format(str(diseases[i])))
    ax.legend(loc='lower right')

# Son boş alt grafiği gizle
fig.delaxes(axs[2, 1]) #fig.delaxes(axs[2, 1])---_> 3 satir 2 sutunken 

# Alt grafikler arasında boşluk ekle
plt.tight_layout()

## save as pdf
plt.savefig("GO_overlapanalysis_assessment3052024.pdf", format="pdf", bbox_inches="tight")
# Grafiği göster
plt.show()   
