import numpy as np
import scipy as scp
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from sklearn import *
from sklearn.svm import *
from sklearn.preprocessing import label_binarize
import random as rn
from sklearn.metrics import roc_curve, auc


fdir='/glycopeptide2/data/users/tgroth/Programs/TCGA_Analysis/BC_Progression/TCGA_data/'
BC_subtypes=['Luminal_A','Luminal_B','Normal','Basal','Her2_positive']
BC_class_glycogene=pd.read_csv(fdir+'BC_class_glycogene.csv',sep=',',index_col=1)
#Filter out samples that don't fit the PAM50 labels:
BC_class_glycogene=BC_class_glycogene.loc[BC_subtypes,:]
BC_class_glycogene=BC_class_glycogene.drop(columns='Unnamed: 0')

#Binarize the labels of all data:
BC_data_labels=label_binarize(list(BC_class_glycogene.index),classes=BC_subtypes)

#Load the PAM50 data
BC_class_allgene=pd.read_csv(fdir+'BC_class_allgene.csv',sep=',',index_col=0)
BC_label_pam50=pd.read_csv(fdir+'BC_label_pam50.tsv',sep='\t',index_col=0)
BC_subtype_pam50=BC_class_allgene[list(BC_label_pam50.columns)]
BC_subtype_pam50=BC_subtype_pam50.loc[BC_subtypes,:]

C_range=np.logspace(-2,2,4)

#Fit a linear SVM for the different Cs, then pick the best-performing one
def Fit_Linear_SVM(C,exp_data,k_fold):
    #Loop through losses, store accuracy statistics:
    acc_data=pd.DataFrame(index=range(len(C)),columns=['C','mean_acc','std_dev_acc'])
    for Cind in range(len(C)):
        Cval=C[Cind]
        print('Cval = %s'%(Cval))
        #Define the classifier:
        clf=LinearSVC(C=Cval)
        #Loop through k_fold iterations,storing accuracy info
        acc_list=[]
        for i in range(k_fold):
            print('Cross %s'%(i))
            n_samples=exp_data.shape[0]
            train_ind=rn.sample(range(n_samples),int(np.round(0.70*n_samples)))
            test_ind=np.setdiff1d(list(range(n_samples)),train_ind)

            train_data=exp_data.iloc[train_ind,:]
            test_data=exp_data.iloc[test_ind,:]

            train_label=list(exp_data.index[train_ind])
            test_label=list(exp_data.index[test_ind])

            #Fit Classifier
            clf.fit(train_data.values,train_label)
            #Make Predictions:
            predict_labels=clf.predict(test_data.values)

            #Find the accuracy
            acc=sum([predict_labels[x]==test_label[x] for x in range(len(test_label))])/len(test_label)
            acc_list.append(acc)
        acc_data.loc[Cind,'C']=Cval
        acc_data.loc[Cind,'mean_acc']=np.mean(acc_list)
        acc_data.loc[Cind,'std_dev_acc']=np.std(acc_list)

    return acc_data
            
def Optim_Fit(C,exp_data,label_data,k_fold):
    clf=SVC(kernel='linear',C=C,max_iter=10000)
    acc_list=[]
    for i in range(k_fold):
        print('Cross %s'%(i))
        n_samples=exp_data.shape[0]
        train_ind=rn.sample(range(n_samples),int(np.round(0.70*n_samples)))
        test_ind=np.setdiff1d(list(range(n_samples)),train_ind)

        train_data=exp_data.iloc[train_ind,:]
        test_data=exp_data.iloc[test_ind,:]

        train_label=list(exp_data.index[train_ind])
        test_label=list(exp_data.index[test_ind])

        #Fit Classifier
        clf.fit(train_data.values,train_label)
        #Make Predictions:
        predict_labels=clf.predict(test_data.values)

        #Find the accuracy
        acc=sum([predict_labels[x]==test_label[x] for x in range(len(test_label))])/len(test_label)
        acc_list.append(acc)

    #Just display the results
    print('Accuracy: %s +/- %s' %(np.mean(acc_list),np.std(acc_list)))

    return clf

def plot_coefficients(clf, feature_names,label_names,top_features=20):
    coef=clf.coef_
    pair_list=[]
    for i in range(len(label_names)-1):
        for j in range(i+1,len(label_names)):
            entry=[label_names[i],label_names[j]]
            pair_list.append(entry)

    for elt in range(len(coef)):
        coef_part=coef[elt]
        top_pos_coef=np.argsort(coef_part)[-top_features:]
        top_neg_coef=np.argsort(coef_part)[:top_features]
        top_coefs=np.hstack([top_neg_coef,top_pos_coef])
        plt.figure()
        colors=['red' if c >0 else 'blue' for c in coef_part[top_coefs]]
        plt.bar(np.arange(2*top_features),coef_part[top_coefs],color=colors)
        feature_names=np.array(feature_names)
        plt.xticks(np.arange(0,(2*top_features)),feature_names[top_coefs],rotation=45,ha='right')
        plt.title('%s vs %s Classification, Significant Features'%(pair_list[elt][0],pair_list[elt][1]))
        plt.show()


def Linear_SVM_ROC(Cval,exp_data,label_names,k_fold):
    #Keeps track of true positive rate / false positive rate:

    #Classifier
    clf=SVC(C=Cval,kernel='linear',probability=True,max_iter=10000)
    
    tprs = []
    fprs = []
    aucs=[]
    mean_fpr = np.linspace(0,1,100)
    i=0

    for cv in range(k_fold):
        print('Cross %s'%(cv))
        n_samples=exp_data.shape[0]
        train_ind=rn.sample(range(n_samples),int(np.round(0.70*n_samples)))
        test_ind=np.setdiff1d(list(range(n_samples)),train_ind)

        train_data=exp_data.iloc[train_ind,:]
        test_data=exp_data.iloc[test_ind,:]

        train_label=list(exp_data.index[train_ind])
        test_label=list(exp_data.index[test_ind])

        #Fit Classifier
        clf.fit(train_data.values,train_label)
        #Make Predictions:
        predict_labels=clf.predict_proba(test_data.values)

        
        fpr = dict()
        tpr = dict()
        tpr_interp=dict()
        roc_auc = dict()
        label_data=label_binarize(test_label,classes=label_names)
        for i in range(len(label_names)):
            fpr[i],tpr[i], _ = roc_curve(label_data[:,i],predict_labels[:,i])
            tpr_interp[i]=np.interp(mean_fpr,fpr[i],tpr[i])
            roc_auc[i]=auc(fpr[i],tpr[i])

        tprs.append(tpr_interp)
        fprs.append(fpr)
        aucs.append(roc_auc)


    mean_tprs=[]
    tpr_uppers=[]
    tpr_lowers=[]
    mean_aucs=[]
    std_aucs=[]
    std_tprs=[]
    for i in range(len(label_names)):
        tpr_class_dct={}
        for x in range(len(tprs)):
            tpr_class_dct[x]=pd.Series(tprs[x][i])
        tpr_class_DF=pd.DataFrame(tpr_class_dct)
        aucs_class_dct={}
        for x in range(len(aucs)):
            aucs_class_dct[x]=pd.Series(aucs[x][i])
            
        aucs_class_DF=pd.DataFrame(aucs_class_dct)
        mean_tprs.append(tpr_class_DF.mean(axis=1).values)
        mean_aucs.append(auc(mean_fpr,mean_tprs[i]))
        std_aucs.append(aucs_class_DF.std(axis=1).values[0])
        std_tprs.append(tpr_class_DF.std(axis=1).values)
        tpr_uppers.append(np.minimum(mean_tprs[i] + std_tprs[i],1))
        tpr_lowers.append(np.maximum(mean_tprs[i] - std_tprs[i],0))

    return mean_tprs,tpr_uppers,tpr_lowers,mean_aucs,std_aucs,std_tprs

# Fit data for PAM50 and Glycogene Groups:


def draw_roc(label_names,plot_index,mean_fpr,mean_tprs,mean_aucs,tpr_lowers,tpr_uppers,std_aucs,std_tprs,
                 ref_mean_tprs,ref_mean_aucs,ref_tpr_lowers,ref_tpr_uppers,ref_std_aucs,ref_std_tprs):
    plt.figure()
    plt.plot(mean_fpr,ref_mean_tprs[plot_index],color='magenta',lw=2,label='PAM50 ROC curve (area = %0.2f)' % ref_mean_aucs[plot_index])
    plt.plot(mean_fpr,mean_tprs[plot_index],color='darkorange',lw=2,label='Glycogene ROC curve (area = %0.2f)' %(mean_aucs[plot_index]))
    plt.fill_between(mean_fpr,ref_tpr_lowers[plot_index],ref_tpr_uppers[plot_index],color='magenta',alpha=0.2,label=r'$\pm$ 1 std.dev')
    plt.fill_between(mean_fpr,tpr_lowers[plot_index],tpr_uppers[plot_index],color='darkorange',alpha=0.2,label=r'$\pm$ 1 std.dev')
    plt.plot([0,1],[0,1],color='navy',lw=2,linestyle='--')
    plt.xlim([-0.1,1.0])
    plt.ylim([0.0,1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positivie Rate')
    plt.title('Breast Cancer Glycogenes to Classify %s (With PAM50 Reference)' % (label_names[plot_index]))
    plt.legend(loc='lower right')
    plt.show()
    
##for elt in range(n_classes):
##    draw_roc(elt,mean_fpr,mean_tprs,mean_aucs,tpr_lowers,tpr_uppers,std_aucs,std_tprs,
##             ref_mean_tprs,ref_mean_aucs,ref_tpr_lowers,ref_tpr_uppers,ref_std_aucs,ref_std_tprs,
##             centroids)
