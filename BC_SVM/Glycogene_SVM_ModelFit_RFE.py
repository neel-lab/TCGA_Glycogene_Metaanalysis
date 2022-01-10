import numpy as np
import scipy as scp
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from sklearn import *
from sklearn.svm import *
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_validation
from sklearn.feature_selection import RFE
import random as rn
from sklearn.metrics import roc_curve, auc
import pickle as pkl


fdir='/glycopeptide2/data/users/tgroth/Programs/TCGA_Analysis/BC_Progression/TCGA_data/'
BC_subtypes=['Luminal_A','Luminal_B','Normal','Basal','Her2_positive']
BC_subtype_integer={
    'Luminal_A':0,
    'Luminal_B':1,
    'Normal':2,
    'Basal':3,
    'Her2_positive':4
    }
BC_class_glycogene=pd.read_csv(fdir+'BC_class_glycogene.csv',sep=',',index_col=1)
#Filter out samples that don't fit the PAM50 labels:
BC_class_glycogene=BC_class_glycogene.loc[BC_subtypes,:]
BC_class_glycogene=BC_class_glycogene.drop(columns='Unnamed: 0')


#Do a grid search on some SVM parameters
# Radial Basis Function (Gaussian) and Linear Models tested:
param_grid=[
    {'C':list(np.logspace(-4,3,20)),'kernel':['linear']},
    ]

svc=svm.SVC()
clf=GridSearchCV(svc,param_grid,cv=10)
label_int=np.array([int(BC_subtype_integer[x]) for x in list(BC_class_glycogene.index)])
clf.fit(BC_class_glycogene.values,label_int)

#Get the best classifier from clf

best_clf = clf.best_estimator_

#Perform multiple recursive feature eliminations on the best model to find consensus gene list
print('Perform Recursive Feature Elimination on the Best Classifier')
rfe=RFE(estimator=best_clf,n_features_to_select=50,step=2)

sig_GG_tally=pd.DataFrame(index=list(BC_class_glycogene.columns),columns=['Frequency'],data=[0]*len(BC_class_glycogene.columns))
RFE_cvNum=1000
for i in range(RFE_cvNum):
    #print('Fit %s' %(i+1))
    rand_samp=rn.sample(range(BC_class_glycogene.shape[0]),int(0.75*(BC_class_glycogene.shape[0])))
    rfe.fit(BC_class_glycogene.iloc[rand_samp,:],label_int[rand_samp])
        
    #Get the feature rankings:
    gg_ranks=rfe.ranking_
    sig_features=[gg_ranks[x]==1 for x in range(len(gg_ranks))]
    sig_GlycoGenes=list(BC_class_glycogene.columns[sig_features])
    for g in sig_GlycoGenes:
        sig_GG_tally.loc[g,]=sig_GG_tally.loc[g,]+1

#sig_GG_tally=sig_GG_tally[sig_GG_tally['Frequency']>0]

#sig_GG_tally.to_csv(fdir+'Bootstrapped_GlycoSVM_list.csv')

#print('All Done!')

#Load the 1000 bootstraps here to avoid running again...:

with open(fdir+'Glycogene_RFE_1000.pkl','rb') as f:
    sig_GG_tally=pkl.load(f)

#Do glycogene titration, see how many glycogenes it takes for accuracy to level off.

#Sort the sig_GG_tally series:
sig_GG_tally=sig_GG_tally.sort_values(axis=0,ascending=False)
#Start with the top 5 most frequently occurring glycogenes in the RFE step.
#  find out how many genes it takes for the accuracy to level off:
import sklearn as sk
titration_means=[]
titration_sds=[]
for i in range(5,sig_GG_tally.shape[0]):
    print('Titrating top '+ str(i) +' genes')
    #Select Genes:
    genes=sig_GG_tally.index[0:i]
    #Extract the glycogene information for training:
    geneData=BC_class_glycogene[genes]

    #Train and test of optimal svm:
    cv_svm_results=sk.model_selection.cross_validate(best_clf,geneData,label_int,cv=10)
    #get the mean and sd of accuracy
    acc_mean=np.mean(cv_svm_results['test_score'])
    acc_sd=np.std(cv_svm_results['test_score'])

    titration_means.append(acc_mean)
    titration_sds.append(acc_sd)


fig = plt.figure()
ax1=fig.add_subplot(111)
ax1.plot(range(5,sig_GG_tally.shape[0]),titration_means,color='red',label='Mean Accuracy')
ax1.set_ylabel('Mean Accuracy')

#Compute linear moving average of the differences: steps of 5
wsize=20
rolling_ave=pd.Series(np.diff(titration_means)).rolling(window=wsize).mean()

ax2=ax1.twinx()
ax2.plot(range(6,sig_GG_tally.shape[0]),rolling_ave,color='green',label=r'$\Delta$ Accuracy')
ax2.set_ylabel(r'$\Delta$ Accuracy Moving Ave')
fig.legend()
fig.show()

# Mean accuracy +/- sd dev
plt.plot(range(5,sig_GG_tally.shape[0]),titration_means,color='red',label='Mean Accuracy')
plt.fill_between(range(5,sig_GG_tally.shape[0]),\
    np.array(titration_means)-np.array(titration_sds),\
                 np.array(titration_means)+np.array(titration_sds),\
                 color='red',alpha=0.2,label=r'$\pm$ 1 std.dev')

plt.title('Glycogene Feature Titration:\nMean Accuracy vs Number of Top Discriminating Glycogenes')
plt.xlabel('Number of Top Glycogenes')
plt.ylabel('Mean Accuracy')
plt.legend(loc='lower right')
plt.axvline(x=30,color='black','--')

#30 looks like a good number of genes to take.  The accuracy levels off after adding this many
# genes.

# Save the top 30 in the list:
sig_GG_tally_best=sig_GG_tally.iloc[0:30,:]

sig_GG_tally_best.to_csv(fdir+'Bootstrapped_GlycoSVM_list_1000.csv')
