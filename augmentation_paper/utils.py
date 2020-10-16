import pandas as pd
import numpy as np
import seaborn as sns
from collections import OrderedDict
import heapq
from Queue import PriorityQueue
from tqdm import tqdm
import pylab,random,cPickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from copy import copy
import time
idx = pd.IndexSlice
from sklearn.metrics import r2_score,mean_squared_error,roc_auc_score, accuracy_score,roc_curve,auc, precision_recall_curve,average_precision_score
from scipy.stats import ranksums,rankdata,pearsonr,spearmanr,mannwhitneyu
from os.path import dirname, basename,join,exists
from os import makedirs,system,listdir,rmdir
from matplotlib.colors import ListedColormap
import multiprocessing as mp
#import dill
import cPickle
import glob
class PopulationCoverage(object):
    def __init__(self,input_epitope,allele_map,country_list):
        self.input_epitope_affinity=input_epitope
        self.allele_map=allele_map
        self.set_country(country_list)
    
    def beam_search(self,beam_size=20,cutoff=0.9,max_round=1000):
        current_max_coverage=0.0
        beams=[]
        curr_beam={}
        cnt=0
        while current_max_coverage<cutoff and cnt<max_round:
            print 'beamsearch round ',cnt
            next_beam={}
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x.split('_') for x in curr_beam.keys()]
            pbar = tqdm(total=len(curr_candidates)*len(self.input_epitope_affinity.index))
            for curr_epitopes in curr_candidates:
                for next_epitope in self.input_epitope_affinity.index:
                    pbar.update(1)
                    if next_epitope in curr_epitopes:
                        continue
                    test_set=sorted(curr_epitopes+[next_epitope])
                    key='_'.join(test_set)
                    #if not key in next_beam:
                    next_beam[key]=self.overall_coverage(test_set,verbose=False)
            pbar.close()
            print 'total combinations for next beam:',len(next_beam)
            curr_beam=OrderedDict(sorted(next_beam.items(), key=lambda t: t[1],reverse = True)[:beam_size])
            current_max_coverage=curr_beam.items()[0][1]
            beams.append(curr_beam)
            print 'current beam: ',curr_beam.items()
            cnt+=1
        print 'Coverage cutoff reached, final solution:',curr_beam.items()[0]
        print 'Per region details:'
        details=self.overall_coverage(epitopes=curr_beam.items()[0][0].split('_'))
        return curr_beam,details,beams

    def overall_coverage(self,epitopes,verbose=True):
        country_coverage=[]
        country_detail={}
        if isinstance(epitopes,str):
            epitopes=[epitopes]
        if len(epitopes)==1:
            binding=self.input_epitope_affinity.loc[epitopes[0]].fillna(0.0)
        else:
            single_binding=self.input_epitope_affinity.loc[epitopes].fillna(0.0)
            binding=1-(1-single_binding).product(axis=0)
        for country in self.country_list:
            #if verbose:
                #print country
            result=self.compute_coverage(country,binding,verbose=verbose)
            country_coverage.append(result[0])
            if verbose:
                #print 'overall coverage:',country_coverage[-1]
                country_detail[country]=[result[0]]+result[1]
        if verbose:
            return np.mean(country_coverage),country_detail
        else:
            return np.mean(country_coverage)
    
    def compute_coverage(self,country,binding,verbose=True):
        locus_binding=[]
        for locus in binding.index.levels[0]:
            single=binding[locus]
            double_binding=1-np.outer(1-single.values,1-single.values)
            np.fill_diagonal(double_binding,single.values)
            locus_binding.append(np.sum(double_binding*(self.diploid_map[country][locus].values)))
        if verbose:
            #print 'locus binding probability:',locus_binding
            return 1.0-np.prod(1-np.asarray(locus_binding)),locus_binding
        return [1.0-np.prod(1-np.asarray(locus_binding))]
        
    def precompute_diploid(self,country_list):
        self.diploid_map={}
        for country in country_list:
            self.diploid_map[country]={}
            for locus in self.allele_map.columns.levels[0]:
                single=self.allele_map.loc[country][locus]
                diploid=pd.DataFrame(np.outer(single.values,single.values),index=single.index,columns=single.index)
                self.diploid_map[country][locus]=diploid

    def set_country(self,country_list):
        self.country_list=country_list
        self.precompute_diploid(country_list)

class KthLargest2(object):
    def __init__(self,k,initial=None):
        self.k=k
        self.heap=[]
        if initial:
            self.heap=initial
        heapq.heapify(self.heap)
        while len(self.heap)>k:
            heapq.heappop(self.heap)
        self.keys=self.get_keys()

    def get_keys(self):
        return [t[1] for t in self.heap]

    def add(self,val):
        if not val[1] in self.keys:
            if len(self.heap)<self.k:
                heapq.heappush(self.heap,val)
            else:
                heapq.heappushpop(self.heap,val)
            self.keys=self.get_keys()
    def return_dict(self):
        return OrderedDict([(k[0],[v,k[1]]) for v,k in sorted(self.heap,reverse=True)])
    def return_dict2(self):
        return [(v,k) for v,k in sorted(self.heap,reverse=True)]
    
class PopulationCoverage2(object):
    def __init__(self,input_epitope,hap_map,country_list,candidates,Sp,pre_map=None,base_counting=None):
        self.input_epitope_affinity=input_epitope #this is a dictionary
        self.hap_map=hap_map
        self.base_dir=base_counting
        self.candidates=candidates
        self.Sp=Sp
        #print base_counting
        if pre_map:
            #self.pre_map=pre_map
            self.set_country(country_list,precompute=False)
        else:
            self.set_country(country_list,precompute=False)

    def overall_multi(self,epitopes_long,lower=0,verbose=False,pre_map=None,return_pre=False,base=None):
        if base:
            self.base_dir=base
            basename=base.split('.pkl')[0].split('_')[-2]
        #print base,basename
        focus=self.candidates.loc[epitopes_long]
        obj=0.0
        result={}
        pept={}
        for mhc in ['MHC1','MHC2']:
            epitopes=set([a for x in focus['compressed_'+mhc].values for a in x.split('_') if len(a)>0])
            pept[mhc]=list(epitopes)
            if base:
                epitopes=epitopes-set(self.Sp[basename+'_'+mhc])
            coverage=self.overall_coverage(epitopes=list(epitopes),lower=lower,pre_map=pre_map,verbose=verbose,mtype=mhc)
            if verbose:
                #print mhc,coverage[0]
                obj+=coverage[0]
                result[mhc]=coverage
            else:
                #print mhc,coverage
                obj+=coverage
                result[mhc]=coverage
        return obj/2.0,result,pept

    def overall_coverage(self,epitopes,lower=0,verbose=True,pre_map=None,return_pre=False,mtype='MHC1',base=None):
        if base:
            self.base_dir=base
            basename=base.split('.pkl')[0].split('_')[-2]
        country_coverage=[]
        country_detail={}
        if isinstance(epitopes,str):
            epitopes=[epitopes]
        counting=self.input_epitope_affinity[mtype].loc[epitopes].fillna(0.0).sum(axis=0)
        for country in self.country_list:
            if pre_map:
                result=self._compute_coverage(country,counting,lower,verbose=verbose,return_pre=return_pre,mtype=mtype)
            else:
                result=self.compute_coverage(country,counting,lower,verbose=verbose,return_pre=return_pre,mtype=mtype)
            country_coverage.append(result[0])
            if verbose:
                country_detail[country]=result
        if verbose:
            return np.mean(country_coverage),country_detail
        else:
            return np.mean(country_coverage)
    
    def compute_coverage(self,country,counting,lower,verbose=True,truncate=True,return_pre=False,savedir=None,mtype='MHC1'):
        global premap
        new_pre3=premap[mtype][country].copy()
        #new_pre3['count']=0
        valid=set(counting[counting>0].index)&self.country_allele[mtype][country]
        prob,hist=self._compute_hist(valid,new_pre3,counting,lower,mtype,return_pre=return_pre,savedir=savedir)
        if verbose:
            if return_pre:
                return prob,hist[0],hist[1]
            else:
                return prob,hist
        else:
            return [prob]

    def _compute_coverage(self,country,counting,lower,verbose=True,truncate=True,return_pre=False,savedir=None,mtype='MHC1'):
        if not self.base_dir is None:
            #print "loading basement counting"
            pre=pd.read_pickle(self.base_dir.format(mtype))
            new_pre3=pre[pre['country']==country]
        else:
            pre=pd.read_pickle('preprocess_haplotype_new{}.pkl'.format(mtype[-1]))
            new_pre3=pre[pre['country']==country]#.drop(labels='country',axis=1)
        #new_pre3['count']=0
        valid=set(counting[counting>0].index)&self.country_allele[mtype][country]
        #print 'searching on %d alleles' % len(valid)
        prob,hist=self._compute_hist(valid,new_pre3,counting,lower,mtype,return_pre=return_pre,savedir=savedir)
        if verbose:
            if return_pre:
                return prob,hist[0],hist[1]
            else:
                return prob,hist
        else:
            return [prob]

    def _compute_hist(self,valid,new_pre3,counting,lower,type,return_pre=False,savedir=None):
        if type=='MHC1':
            for al in valid:
                if al[4]=='A':
                    try:
                        new_pre3.loc[al,'count']=new_pre3.loc[al,'count'].values+counting[al]
                    except:
                        print al
                    try:
                        new_pre3.loc[idx[:,al],'count']=new_pre3.loc[idx[:,al],'count'].values+counting[al]
                    except:
                        pass#print al
                elif al[4]=='B':
                    try:
                        new_pre3.loc[idx[:,:,al],'count']=new_pre3.loc[idx[:,:,al],'count'].values+counting[al]
                    except:
                        pass
                    try:
                        new_pre3.loc[idx[:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
                else:
                    try:
                        new_pre3.loc[idx[:,:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
                    try:
                        new_pre3.loc[idx[:,:,:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
        else:
            for al in valid:
                if 'DRB' in al:
                    try:
                        new_pre3.loc[al,'count']=new_pre3.loc[al,'count'].values+counting[al]
                    except:
                        print al
                    try:
                        new_pre3.loc[idx[:,al],'count']=new_pre3.loc[idx[:,al],'count'].values+counting[al]
                    except:
                        pass#print al
                elif 'HLA-DP' in al:
                    try:
                        new_pre3.loc[idx[:,:,al],'count']=new_pre3.loc[idx[:,:,al],'count'].values+counting[al]
                    except:
                        pass
                    try:
                        new_pre3.loc[idx[:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
                else:
                    try:
                        new_pre3.loc[idx[:,:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
                    try:
                        new_pre3.loc[idx[:,:,:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
        hist=new_pre3.groupby('count').sum().reset_index()
        prob=hist[hist['count']>lower]['freq'].sum()
        if return_pre:
            return prob,(hist,new_pre3)
        else:
            return prob,hist

    def set_country(self,country_list,precompute=True):
        self.country_list=country_list
        self.country_allele={'MHC1':{},'MHC2':{}}
        if precompute:
            self.precompute_hap(country_list)
        for country in country_list:
            for mhc in self.hap_map:
                v=self.hap_map[mhc].loc[country][self.hap_map[mhc].loc[country]>0].index.values
                self.country_allele[mhc][country]=set([x for it in v for x in it])

import pathos.pools as pp
class KthLargest(object):
    def __init__(self,k,initial=None):
        self.k=k
        self.heap=[]
        if initial:
            self.heap=initial
        heapq.heapify(self.heap)
        while len(self.heap)>k:
            heapq.heappop(self.heap)
        self.keys=self.get_keys()

    def get_keys(self):
        return [t[1] for t in self.heap]

    def add(self,val):
        if not val[1] in self.keys:
            if len(self.heap)<self.k:
                heapq.heappush(self.heap,val)
            else:
                heapq.heappushpop(self.heap,val)
            self.keys=self.get_keys()
    def return_dict(self):
        return OrderedDict([(k,v) for v,k in sorted(self.heap,reverse=True)])
        
class PopulationCoverage3(object):
    def __init__(self,input_epitope,hap_map,country_list,outdir='test',pre_map=None,base_counting=None):
        self.input_epitope_affinity=input_epitope
        self.hap_map=hap_map
        self.outdir=outdir
        self.base_dir=base_counting
        print base_counting
        if pre_map:
            #self.pre_map=pre_map
            self.set_country(country_list,precompute=False)
        else:
            self.set_country(country_list,precompute=False)

    def overall_coverage(self,epitopes,lower=0,verbose=True,pre_map=None,return_pre=False,typem='mhc1_haplotype',base=None):
        if base:
            self.base_dir=base
        country_coverage=[]
        country_detail={}
        if isinstance(epitopes,str):
            epitopes=[epitopes]
        counting=self.input_epitope_affinity.loc[epitopes].fillna(0.0).sum(axis=0)
        for country in self.country_list:
            if pre_map:
                result=self._compute_coverage(country,counting,lower,verbose=verbose,return_pre=return_pre,typem=typem)
            else:
                result=self.compute_coverage(country,counting,lower,verbose=verbose,return_pre=return_pre,typem=typem)
            country_coverage.append(result[0])
            if verbose:
                country_detail[country]=result
        if verbose:
            return np.mean(country_coverage),country_detail
        else:
            return np.mean(country_coverage)
    
    def compute_coverage(self,country,counting,lower,verbose=True,truncate=True,return_pre=False,savedir=None,typem='mhc1_haplotype'):
        global premap
        new_pre3=premap[country].copy()
        #new_pre3['count']=0
        valid=set(counting[counting>0].index)&self.country_allele[country]
        prob,hist=self._compute_hist(valid,new_pre3,counting,lower,typem,return_pre=return_pre,savedir=savedir)
        if verbose:
            if return_pre:
                return prob,hist[0],hist[1]
            else:
                return prob,hist
        else:
            return [prob]

    def _compute_coverage(self,country,counting,lower,verbose=True,truncate=True,return_pre=False,savedir=None,typem='mhc1_haplotype'):
        if not self.base_dir is None:
            print "loading basement counting"
            pre=pd.read_pickle(self.base_dir)
            new_pre3=pre[pre['country']==country]
        else:
            pre=pd.read_pickle('preprocess_haplotype_new{}.pkl'.format(typem.split('_')[0][-1]))
            new_pre3=pre[pre['country']==country]#.drop(labels='country',axis=1)
        #new_pre3['count']=0
        valid=set(counting[counting>0].index)&self.country_allele[country]
        #print 'searching on %d alleles' % len(valid)
        prob,hist=self._compute_hist(valid,new_pre3,counting,lower,typem,return_pre=return_pre,savedir=savedir)
        if verbose:
            if return_pre:
                return prob,hist[0],hist[1]
            else:
                return prob,hist
        else:
            return [prob]

    def _compute_hist(self,valid,new_pre3,counting,lower,type,return_pre=False,savedir=None):
        if type=='mhc1_haplotype':
            for al in valid:
                if al[4]=='A':
                    try:
                        new_pre3.loc[al,'count']=new_pre3.loc[al,'count'].values+counting[al]
                    except:
                        print al
                    try:
                        new_pre3.loc[idx[:,al],'count']=new_pre3.loc[idx[:,al],'count'].values+counting[al]
                    except:
                        pass#print al
                elif al[4]=='B':
                    try:
                        new_pre3.loc[idx[:,:,al],'count']=new_pre3.loc[idx[:,:,al],'count'].values+counting[al]
                    except:
                        pass
                    try:
                        new_pre3.loc[idx[:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
                else:
                    try:
                        new_pre3.loc[idx[:,:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
                    try:
                        new_pre3.loc[idx[:,:,:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
        else:
            for al in valid:
                if 'DRB' in al:
                    try:
                        new_pre3.loc[al,'count']=new_pre3.loc[al,'count'].values+counting[al]
                    except:
                        print al
                    try:
                        new_pre3.loc[idx[:,al],'count']=new_pre3.loc[idx[:,al],'count'].values+counting[al]
                    except:
                        pass#print al
                elif 'HLA-DP' in al:
                    try:
                        new_pre3.loc[idx[:,:,al],'count']=new_pre3.loc[idx[:,:,al],'count'].values+counting[al]
                    except:
                        pass
                    try:
                        new_pre3.loc[idx[:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
                else:
                    try:
                        new_pre3.loc[idx[:,:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
                    try:
                        new_pre3.loc[idx[:,:,:,:,:,al],'count']=new_pre3.loc[idx[:,:,:,:,:,al],'count'].values+counting[al]
                    except:
                        pass
        hist=new_pre3.groupby('count').sum().reset_index()
        prob=hist[hist['count']>lower]['freq'].sum()
        if return_pre:
            return prob,(hist,new_pre3)
        else:
            return prob,hist
        
    def precompute_hap(self,country_list):
        self.pre_map={}
        for country in country_list:
            self.pre_map[country]={'alleles':[],'freq':[]}
            single=self.hap_map.loc[country]
            test=single[single>0].index.values
            print 'precomputing haplotype combination for '+country
            with tqdm(total=(len(test)*(len(test)-1)/2+len(test))) as pbar:
                for i,hap1 in enumerate(test):
                    for j in range(i,len(test)):
                        hap2=test[j]
                        pbar.update(1)
                        self.pre_map[country]['freq'].append(hap[hap1].loc[country]*hap[hap2].loc[country]*(2-(i==j)))
                        self.pre_map[country]['alleles'].append(np.union1d(hap1,hap2))
            self.pre_map[country]=pd.DataFrame(self.pre_map[country])

    def set_country(self,country_list,precompute=True):
        self.country_list=country_list
        self.country_allele={}
        if precompute:
            self.precompute_hap(country_list)
        for country in country_list:
            v=self.hap_map.loc[country][self.hap_map.loc[country]>0].index.values
            self.country_allele[country]=set([x for it in v for x in it])
