import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import heapq
from collections import OrderedDict
from os.path import dirname, basename,join,exists
from os import makedirs,system,listdir,rename
from tqdm import tqdm
from argparse import ArgumentParser
import sys,glob
import cPickle
from matplotlib.colors import ListedColormap
from datetime import datetime
import pathos.pools as pp
import time
import random
from copy import copy
idx = pd.IndexSlice
import warnings
warnings.filterwarnings("ignore")

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
        return OrderedDict([(k[0],[v,k[1]]) for v,k in sorted(self.heap,reverse=True)])
    def return_dict2(self):
        return [(v,k) for v,k in sorted(self.heap,reverse=True)]

def diff1d(curr,seq,cut=3):
    for x in curr:
        if (x in seq) or (seq in x):
            if abs(len(x)-len(seq))<=cut:
                return False
        else:
            for diff_l in range(1,cut):
                for diff_r in range(1,cut-diff_l+1):
                    if x[diff_l:]==seq[:-diff_r] or seq[diff_l:]==x[:-diff_r]:
                        return False
    return True

def initializer():
    global premap
    pre=pd.read_pickle('preprocess_haplotype_new1.pkl')
    premap1={}
    for c in ['White','Black','Asians']:
        premap1[c]=pre[pre['country']==c].drop(labels='country',axis=1)
    pre=pd.read_pickle('preprocess_haplotype_new2.pkl')
    premap2={}
    for c in ['White','Black','Asians']:
        premap2[c]=pre[pre['country']==c].drop(labels='country',axis=1)
    premap={'mhc1':premap1,'mhc2':premap2}

class PopulationCoverage(object):
    def __init__(self,input_epitope,hap_map,country_list,outdir,pre_map=None,base_counting=None):
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

    def overall_coverage(self,epitopes,lower=0,verbose=True,pre_map=None,return_pre=False,mtype='mhc1'):
        country_coverage=[]
        country_detail={}
        if isinstance(epitopes,str):
            epitopes=[epitopes]
        counting=self.input_epitope_affinity.loc[epitopes].fillna(0.0).sum(axis=0)
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
    
    def compute_coverage(self,country,counting,lower,verbose=True,truncate=True,return_pre=False,savedir=None,mtype='mhc1'):
        global premap
        new_pre3=premap[country].copy()
        #new_pre3['count']=0
        valid=set(counting[counting>0].index)&self.country_allele[country]
        prob,hist=self._compute_hist(valid,new_pre3,counting,lower,mtype,return_pre=return_pre,savedir=savedir)
        if verbose:
            if return_pre:
                return prob,hist[0],hist[1]
            else:
                return prob,hist
        else:
            return [prob]

    def _compute_coverage(self,country,counting,lower,verbose=True,truncate=True,return_pre=False,savedir=None,mtype='mhc1'):
        if not self.base_dir is None:
            print "loading basement counting"
            pre=pd.read_pickle(self.base_dir)
            new_pre3=pre[pre['country']==country]
        else:
            pre=pd.read_pickle('preprocess_haplotype_new{}.pkl'.format(mtype[-1]))
            new_pre3=pre[pre['country']==country]#.drop(labels='country',axis=1)
        #new_pre3['count']=0
        valid=set(counting[counting>0].index)&self.country_allele[country]
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
        if type=='mhc1':
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
        self.country_allele={}
        if precompute:
            self.precompute_hap(country_list)
        for country in country_list:
            v=self.hap_map.loc[country][self.hap_map.loc[country]>0].index.values
            self.country_allele[country]=set([x for it in v for x in it])
    
class PopulationCoverage2(object):
    def __init__(self,input_epitope,hap_map,country_list,outdir,candidates,pre_map=None,base_counting=None):
        self.input_epitope_affinity=input_epitope #this is a dictionary
        self.hap_map=hap_map
        self.outdir=outdir
        self.base_dir=base_counting
        self.candidates=candidates
        print base_counting
        if pre_map:
            #self.pre_map=pre_map
            self.set_country(country_list,precompute=False)
        else:
            self.set_country(country_list,precompute=False)

    def beam_search_parallel(self,beam_size=20,cutoff=0.9,min_cutoff=5,\
                            max_round=20,curr_beam={},curr_min=0,nworker=20,bs=50,diverse_cut=3,augment=False,upper=None,lower=None):
        print 'Using %d workers' % nworker
        if len(curr_beam)>0:
            current_max_coverage,details=self.overall_multi(epitopes_long=curr_beam.items()[0][0].split('_'),lower=curr_min,pre_map=True,verbose=False)
        else:
            current_max_coverage,details=self.overall_multi(epitopes_long=[],lower=curr_min,pre_map=True,verbose=False)
        beams=[]
        curr_length=max([len(k[0].split('_')) for k in curr_beam.items()]) if len(curr_beam)>0 else 0
        #curr_length=len(curr_beam.items()[0][0].split('_')) if len(curr_beam)>0 else 0
        outdir=join(self.outdir,'plots')
        print 'current beam: ',curr_beam.items()
        print 'current lower bound: ',curr_min
        if args.initial_cut:
            iter_cutoff=args.initial_cut
        else:
            if augment:
                if curr_min in range(0,min_cutoff,4)+[min_cutoff]:
                    iter_cutoff=max(upper['average'].iloc[curr_min]-0.01,(args.ratio*upper['average'].iloc[curr_min]+(1-args.ratio)*current_max_coverage))
                else:
                    iter_cutoff=max(lower['average'].iloc[curr_min]+0.005,(args.ratio_low*max(current_max_coverage,lower['average'].iloc[curr_min])+(1-args.ratio_low)*upper['average'].iloc[curr_min]))
            else:
                if curr_min>0:
                    if curr_min in range(0,min_cutoff,4)+[min_cutoff]:
                        iter_cutoff=max(upper['average'].iloc[curr_min]-0.01,(args.ratio*upper['average'].iloc[curr_min]+(1-args.ratio)*current_max_coverage))
                    else:
                        iter_cutoff=args.ratio_low*current_max_coverage+(1-args.ratio_low)*upper['average'].iloc[curr_min]
                else:
                    if args.type.split('_')[0][-1]=='1':
                        iter_cutoff=max(0.99,cutoff)
                    else:
                        iter_cutoff=max(0.99,cutoff)
        old_max=current_max_coverage
        while (current_max_coverage<iter_cutoff or curr_min<min_cutoff) and curr_length<max_round:
            pool=pp._ProcessPool(processes=nworker,initializer=initializer,initargs=())
            print 'beamsearch round ',curr_length,'current min:', curr_min, 'cutoff', iter_cutoff,'curr_max',current_max_coverage
            t0 = time.time()
            self.next_beam=KthLargest(k=beam_size)
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x.split('_') for x in curr_beam.keys()]
            all_args=[]
            #print 'current candidates:',curr_candidates
            test_epitopes=[]
            test_epitopes_backup=[]
            seen_curr={}
            test_epitopes=[x+[y] for x in curr_candidates for y in self.candidates.index if not y in x]
            for divide in range(10):
                bs=len(test_epitopes)//((divide+1)*nworker)+int(len(test_epitopes)%((divide+1)*nworker)>0)
                if bs<80:
                    print 'using batch size =',bs
                    break
            #test_epitopes=[x+[y] for x in curr_candidates for y in self.input_epitope_affinity.index if diff1d(x,y,cut=diverse_cut)]
            for i in range(0,len(test_epitopes),bs):
                all_args.append((copy(test_epitopes[i:min(i+bs,len(test_epitopes))]),curr_min,beam_size))
            print len(test_epitopes),len(all_args)
            del test_epitopes
            r=pool.map_async(self.batch_scan,all_args,callback=self.batch_update)
            r.get()
            curr_beam=self.next_beam.return_dict()
            current_max_coverage=curr_beam.items()[0][1][0]
            beams.append(curr_beam)
            # check histogram here
            curr_length=max([len(k[0].split('_')) for k in curr_beam.items()])
            #curr_length=len(curr_beam.items()[0][0].split('_'))
            current_median=np.median([v[1][0] for v in curr_beam.items()])
            old_min=curr_min
            if current_median > iter_cutoff or current_max_coverage-old_max<1e-5:
                if current_max_coverage-old_max<1e-5:
                    print "converge, can not futher optimize!"
                else:
                    print 'median min_coverage of {} reached {}, raising min_coverage'.format(curr_min,current_median)
                curr_min+=1
                current_max_coverage,details=self.overall_multi(epitopes_long=curr_beam.items()[0][0].split('_'),lower=curr_min,pre_map=True,verbose=False)
                if augment:
                    if curr_min in range(0,min_cutoff,4)+[min_cutoff]:
                        iter_cutoff=max(upper['average'].iloc[curr_min]-0.01,(args.ratio*upper['average'].iloc[curr_min]+(1-args.ratio)*current_max_coverage))
                    else:
                        iter_cutoff=max(lower['average'].iloc[curr_min]+0.005,(args.ratio_low*max(current_max_coverage,lower['average'].iloc[curr_min])+(1-args.ratio_low)*upper['average'].iloc[curr_min]))
                else:
                    if curr_min>0:
                        if curr_min in range(0,min_cutoff,4)+[min_cutoff]:
                            iter_cutoff=max(upper['average'].iloc[curr_min]-0.01,(args.ratio*upper['average'].iloc[curr_min]+(1-args.ratio)*current_max_coverage))
                        else:
                            iter_cutoff=args.ratio_low*current_max_coverage+(1-args.ratio_low)*upper['average'].iloc[curr_min]
                    else:
                        if args.type.split('_')[0][-1]=='1':
                            iter_cutoff=max(0.99,cutoff)
                        else:
                            iter_cutoff=max(0.99,cutoff)
            old_max=current_max_coverage
            print 'current beam: ',curr_beam.items()
            print 'current lower bound: ',curr_min,iter_cutoff,current_max_coverage
            save_pickle(join(self.outdir,'beam_'+str(curr_length-1)+'.p'),curr_beam.items())
            save_beams(join(self.outdir,'beam_'+str(curr_length-1)),curr_beam,curr_min,old_min)
            t1 = time.time()
            print('time passed: {}'.format(t1-t0))
            pool.close()

        print 'Coverage cutoff reached, final solution:',curr_beam.items()[0]
        print 'Per region details:'
        details=self.overall_coverage(epitopes=curr_beam.items()[0][0].split('_'),lower=min_cutoff,pre_map=True,verbose=True)
        # result_hist=compute_probs(details[1])
        # result_hist.to_csv(join(outdir,'best_result_histogram.csv'))
        np.savetxt(join(outdir,'final_sequences.txt'),curr_beam.items()[0][0].split('_'),fmt='%s')
        return curr_beam.items()[0],details,beams

    def overall_multi(self,epitopes_long,lower=0,verbose=False,pre_map=None,return_pre=False):
        focus=self.candidates.loc[epitopes_long]
        obj=0.0
        result={}
        for mhc in ['mhc1','mhc2']:
            epitopes=set([a for x in focus['compressed_'+mhc].values for a in x.split('_') if len(a)>0])
            try:
                self.input_epitope_affinity[mhc].loc[epitopes]
            except:
                epitopes=[]
            coverage=self.overall_coverage(epitopes=list(epitopes),lower=lower,pre_map=pre_map,verbose=verbose,mtype=mhc)
            if verbose:
                #print mhc,coverage[0]
                obj+=coverage[0]
                result[mhc]=coverage
            else:
                #print mhc,coverage
                obj+=coverage
                result[mhc]=coverage
        return obj/2.0,result


    def batch_scan2(self,inputs):
        global premap
        local_beam=KthLargest(k=5)
        test_sets,lower_bound = inputs
        result=[]
        for test_set in test_sets:
            test_set=sorted(test_set)
            key='_'.join(test_set)
            coverage=self.overall_coverage(test_set,lower_bound,verbose=False)
            result.append((coverage,key))
        for item in result:
                local_beam.add(item)
        return [(x[1],x[0]) for x in local_beam.return_dict().items()]#result
    
    def batch_scan(self,inputs):
        global premap
        test_sets,lower_bound,beamsize = inputs
        local_beam=KthLargest(k=beamsize)
        result=[]
        #print (test_sets[0])
        if not self.base_dir is None:
            #premap=pd.read_pickle(self.base_dir)
            premap={}
            for mhc in ['mhc1','mhc2']:
                pre=pd.read_pickle(self.base_dir.format(mhc))
                premap1={}
                for c in ['White','Black','Asians']:
                    premap1[c]=pre[pre['country']==c].drop(labels='country',axis=1)
                premap[mhc]=premap1
        for test_set in test_sets:
            test_set=sorted(test_set)
            key='_'.join(test_set)
            coverage,detail=self.overall_multi(test_set,lower_bound,verbose=False)
            result.append((coverage,[key,detail]))
        for item in result:
            local_beam.add(item)
        return local_beam.return_dict2()

    def batch_update(self,results):
        print 'called update'
        print len(results[0]),len(results)
        for rs in results:
            for item in rs:
                self.next_beam.add(item)

    def overall_coverage(self,epitopes,lower=0,verbose=True,pre_map=None,return_pre=False,mtype='mhc1'):
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
    
    def compute_coverage(self,country,counting,lower,verbose=True,truncate=True,return_pre=False,savedir=None,mtype='mhc1'):
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

    def _compute_coverage(self,country,counting,lower,verbose=True,truncate=True,return_pre=False,savedir=None,mtype='mhc1'):
        if not self.base_dir is None:
            print "loading basement counting"
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
        if type=='mhc1':
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
        self.country_allele={'mhc1':{},'mhc2':{}}
        if precompute:
            self.precompute_hap(country_list)
        for country in country_list:
            for mhc in self.hap_map:
                v=self.hap_map[mhc].loc[country][self.hap_map[mhc].loc[country]>0].index.values
                self.country_allele[mhc][country]=set([x for it in v for x in it])

def get_parser():
    """Get parser object for script calculate_population_coverage.py."""

    parser = ArgumentParser()

    req_argument = parser.add_argument_group('required arguments')
    parser.add_argument("-m","--method", type=str, default='netmhc',
                        help="which prediction model to use")
    parser.add_argument("-p","--prediction", type=str, default='pred_affinity',
                        help="which type of output metric to use")
    parser.add_argument("-t","--type", type=str, default='mhc1_haplotype',
                        help="which mhc type to use")
    parser.add_argument("-b","--binarize",dest='binary_cutoff', type=float,default=0.638,
                        help="Cutoff for binarizing")
    parser.add_argument("-tr","--truncate",dest='truncate_cutoff', type=float,
                        help="Cutoff for truncating")
    parser.add_argument("-gl","--glyco",dest='glyco_cutoff', type=float,default=0,
                        help="Cutoff for glycolysation")
    parser.add_argument("-mt","--mutation",dest='mutation_cutoff', type=float,default=0.0,
                        help="Cutoff for mutation")
    parser.add_argument("-s","--size",dest='beam_size', type=int, default=1,
                        help="Size of the beam")
    parser.add_argument("-c","--cutoff",dest='coverage_cutoff', type=float, default=0.90,
                        help="Target coverage lower bond when stop beam search")
    parser.add_argument("-oc","--givencut",dest='given_cut', type=float, default=None,
                        help="Target coverage lower bond when stop beam search")
    parser.add_argument("-ic","--initcutoff",dest='initial_cut', type=float,
                        help="Target coverage lower bond when stop beam search")
    parser.add_argument("-high","--ratio_high",dest='ratio', type=float,default=0.85,
                        help="Target coverage upper/lower ratio")
    parser.add_argument("-low","--ratio_low",dest='ratio_low', type=float,default=0.40,
                        help="Target coverage lower/upper ratio")
    parser.add_argument("-lo","--lower",dest='lower_bound', type=int, default=3,
                        help="Target coverage lower bond #peptide when stop beam search")
    parser.add_argument("-r","--maxround",dest='max_round', type=int, default=100,
                        help="Max number of peptides to include")
    parser.add_argument("-f","--frequency",dest='freq_file', type=str, default='iedb',
                        help="which frequency file to use")
    parser.add_argument("-o","--outdir", type=str, default='result',
                        help="path for results.")
    parser.add_argument("-w","--nworker", type=int,
                        help="Number of workers if use parallelization.")
    parser.add_argument("-bs","--batchsize",dest='batchsize', type=int,default=50,
                        help="Number of workers if use parallelization.")
    parser.add_argument("-re","--regions", type=str,
                        default="regions_list_noca.txt",
                        help="filename that contains regions")
    parser.add_argument("-pr","--protein", type=str,
                        help="which proteins the peptides are coming from")
    parser.add_argument("-ba","--basement", type=str,
                        help="which proteins are used as base to augment from")
    parser.add_argument("-bd","--basefile", type=str,
                        help="file name to load basement peptides")
    parser.add_argument("-d","--diversity", type=int, default=3,
                        help="Distance for removing similar sliding windows")
    parser.add_argument("--unroll", action='store_false',
                        help="path for results.")
    parser.add_argument("--downsamp", action='store_true',
                        help="path for results.")
    parser.add_argument("--correction", action='store_true',
                        help="path for results.")
    parser.add_argument("--restart", action='store_true',
                        help="If restarting from previous results or not")
    parser.add_argument("--skippre", action='store_true',
                        help="If skip preprocessing the prediction file or not")
    return parser

def save_detail(fname,result,has_overall=False):
    max_cov=result[0]
    max_detail=result[1]
    with open(fname,'w') as f:
        if has_overall:
            f.writelines('overall coverage:%f\n' %max_cov)
        for country in max_detail:
            f.writelines('%s\t%f\n' % (country,max_detail[country][0]))

def save_beams(fname,beam,curr_min,old_min):
    with open(fname,'w') as f:
        f.writelines('%s\t%d\n' % ('Next lower bound',curr_min))
        f.writelines('%s\t%d\n' % ('Current lower bound',old_min))
        for seq,val in beam.items():
            f.writelines('%s\t%f\t%f\t%f\n' % (seq,val[0],val[1]['mhc1'],val[1]['mhc2']))

def save_pickle(fname,data):
    with open(fname,'w') as f:
        cPickle.dump(data,f)

def load_pickle(fname):
    with open(fname,'rb') as f:
        data=cPickle.load(f)
    return data

def plot_hist(hist_detail,title,savename):
    f,axs=plt.subplots(1,len(hist_detail),figsize=(3.5*len(hist_detail),2.5))
    for j,cntry in enumerate(hist_detail):
        hist=hist_detail[cntry][1]
        axs[j].bar(hist['count'].values,hist['freq'].values,width=1)
        axs[j].set_title(cntry+' '+title)
        axs[j].set_xlabel('# peptide-hla association')
        axs[j].set_ylabel('frequency')
    plt.tight_layout()
    f.savefig(savename, dpi=200,bbox_inches = "tight")

def compute_probs(hist_detail):
    hist_all={'White':[],'Black':[],'Asians':[]}
    for j,cntry in enumerate(hist_detail):
        hist=hist_detail[cntry][1]
        for cut in range(40):
            hist_all[cntry].append(hist[hist['count']>cut]['freq'].sum())
    hist_all=pd.DataFrame(hist_all)
    hist_all['average']=hist_all.mean(axis=1)
    print hist_all
    return hist_all

def plot_detail(model,epitopes,pc,outdir,suffix='',norm=True):
    for epitope in epitopes:
        print epitope
        cv,det=pc.overall_coverage(epitopes=[epitope])
        f, (a0, a1) = plt.subplots(1, 2,figsize=(10,4), gridspec_kw={'width_ratios': [1.5, 2.5]})
        plot_country(det,epitope,a0)
        plot_allele(model,epitope,epitope,a1,norm)
        plt.tight_layout()
        f.savefig(join(outdir,epitope+suffix+'.png'), dpi=200,bbox_inches = "tight")
    f, a0 = plt.subplots(1, 1,figsize=(5.5,4))
    cv,det=pc.overall_coverage(epitopes=epitopes)
    plot_country(det,'Optimal peptide set',a0)
    plt.tight_layout()
    f.savefig(join(outdir,'all_peptide'+suffix+'.png'), dpi=200,bbox_inches = "tight")

def plot_country(country_detail,title,axes):
    d1={'Region':[],'Overall coverage':[]}
    for i,cntry in enumerate(country_detail):
        d1['Overall coverage'].append(country_detail[cntry][0])
        d1['Region'].append(cntry)
    d1=pd.DataFrame(d1)
    d1.loc[i+1]=d1.mean(axis=0)
    d1.loc[i+1,'Region']='Average'
    ax=sns.barplot(x='Region',y='Overall coverage',data=d1,ax=axes)
    ax.set_ylim(0,1)
    for item in ax.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(12)
    axes.set_title(title)

def plot_allele(model,peptide,title,axes,norm=True):
    d=model.loc[peptide].reset_index()
    d.columns=['loci','genotype','Predicted binding']
    ax=sns.barplot(x='genotype',y='Predicted binding',data=d,hue='loci',ax=axes)
    if norm:
        ax.set_ylim(0,1)
    for item in ax.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(6)
    axes.set_title(title)
    axes.legend(ncol=3,loc=0,fancybox=True, framealpha=0.3)

def plot_detail2(model,epitopes,pc,outdir,name='OptiVax(Ours)',suffix='',pal="husl",size=(12,3.5),lgd=True,processed='(truncated)'):
    f, (a0, a1) = plt.subplots(1, 2,figsize=size, gridspec_kw={'width_ratios': [1.5, 2.5]})
    cv,det=pc.overall_coverage(epitopes=epitopes,lower=args.lower_bound,pre_map=True)
    plot_country(det,name+' peptide set',a0)
    #print cv
    model.loc[epitopes].T.fillna(0).plot(kind='bar',stacked=True,\
                               colormap=ListedColormap(sns.color_palette(pal, len(epitopes))),ax=a1)    
    a1.legend(bbox_to_anchor=(1, 1),ncol=len(epitopes)//13+1)
    if not lgd:
        a1.get_legend().remove()
    a1.set_title('Predicted binding'+processed)
    a1.set_xlabel('Allele')
    for item in a1.get_xticklabels():
        item.set_fontsize(8)
    plt.tight_layout()
    f.savefig(join(outdir,'stacked_detail'+suffix+'.png'), dpi=200,bbox_inches = "tight")

def plot_detail3(model,epitopes,pc,outdir,name='OptiVax(Ours)',suffix='',pal="husl",size=(12,3.5),lgd=True,processed='(truncated)'):
    f, (a0, a1) = plt.subplots(1, 2,figsize=size, gridspec_kw={'width_ratios': [1.2,3.2]})
    cv,det=pc.overall_coverage(epitopes=epitopes,lower=args.lower_bound,pre_map=True)
    plot_country(det,name+' peptide set {:.2f}% coverage'.format(cv*100),a0)
    detail=model.loc[epitopes].copy()
    pos=objectives['protein'].loc[epitopes].values
    rename=['{}({})'.format(epitopes[i],pos[i]) for i in range(len(pos))]
    detail.index=rename
    detail.T.fillna(0).plot(kind='bar',stacked=True,\
                               colormap=ListedColormap(sns.color_palette(pal, len(epitopes))),ax=a1)    
    lg=a1.legend(bbox_to_anchor=(1, 1),ncol=len(epitopes)//14+1)
    for text in lg.get_texts():
        if objectives['glyco_probs'].loc[text.get_text().split('(')[0]]>0:
            plt.setp(text, color = 'r')
    if not lgd:
        a1.get_legend().remove()
    a1.set_title('Predicted binding'+processed)
    a1.set_xlabel('Allele')
    for item in a1.get_xticklabels():
        item.set_fontsize(6)
    plt.tight_layout()
    f.savefig(join(outdir,'stacked_detail'+suffix+'.png'), dpi=200,bbox_inches = "tight")
    f2=plt.figure(figsize=(4,3))
    ax=sns.countplot(x='protein',data=objectives.loc[epitopes])
    ax.set_xlabel('Peptide origin',fontsize=12)
    ax.set_ylabel('Count',fontsize=12)
    plt.tight_layout()
    f2.savefig(join(outdir,'protein_distribution'+suffix+'.png'), dpi=200,bbox_inches = "tight")

if __name__ == "__main__":
    args = get_parser().parse_args()
    random.seed(0)
    if args.type=='mhc1':
        if args.freq_file=='iedb':
            frequency=pd.read_pickle('IEDB_population_frequency102_normalized.pkl')
        else:
            print("frequency file type unknown")
            sys.exit(0)
    elif args.type=='mhc2':
        if args.freq_file=='iedb':
            frequency=pd.read_pickle('IEDB_population_frequency_mhc2_normalized.pkl')
        else:
            print("frequency file type unknown")
            sys.exit(0)
    elif 'haplotype' in args.type:
        if 'mhc1' in args.type:
            frequency=pd.read_pickle('haplotype_frequency_marry.pkl')
        elif 'mhc2' in args.type:
            frequency=pd.read_pickle('haplotype_frequency_marry2.pkl')
        elif 'multi' in args.type:
            frequency1=pd.read_pickle('haplotype_frequency_marry.pkl')
            frequency2=pd.read_pickle('haplotype_frequency_marry2.pkl')
            frequency={'mhc1':frequency1,'mhc2':frequency2}
    else:
        print("prediction output type unknown")
        sys.exit(0)

    # if not exists(args.regions):
    #     print("Region file not exist: %s" % args.regions)
    #     sys.exit(0)
    regions=['White','Black','Asians']#np.loadtxt(args.regions,dtype='S',delimiter='\n')

    # if args.method not in ['puffin','deepligand','mhcflurry','netmhc','test','average','mean','netmhcii-3.2','netmhcii-4.0']:
    #     print("prediction method unknown")
    #     sys.exit(0)
    # if args.prediction not in ['pred_prob','ic50','pred_affinity','rank']:
    #     print("prediction output type unknown")
    #     sys.exit(0)

    pred_file1='_'.join(['mhc1_haplotype', 'ensemb-adapt', 'pred_affinity', 'pivot.pkl.gz'])
    pred_file2='_'.join(['mhc2_haplotype', 'netmhcii-4.0-adapt', 'pred_affinity', 'pivot.pkl.gz'])
    if not exists(pred_file1):
        print("prediction file not exist %s" % pred_file)
        sys.exit(0)
    print('reading prediction file: %s,%s' % (pred_file1, pred_file2))
    pred_files={'mhc1':pred_file1,'mhc2':pred_file2}
    
    if not exists(args.outdir):
        makedirs(args.outdir)
        curr_beam={}
        curr_min=0
    else:
        if not args.restart and exists(join(args.outdir,'best_result.txt')):
            print("Result file already exist: %s" % (args.outdir))
            sys.exit(0)
        if args.restart:
            files=glob.glob(join(args.outdir,'beam_*.p'))
            bnum=[int(x.split('_')[-1][:-2]) for x in files]
            if max(bnum)>=args.max_round-1:
                print("Optimization max round is already reached previously")
                if exists(join(args.outdir,'best_result.txt')):
                    if max(bnum)==args.max_round-1:
                        print("Result already exist")
                        sys.exit(0)
                    else:
                        tmstamp=datetime.now().strftime("%Y%m%d%H%M%S")
                        rename(join(args.outdir,'best_result.txt'),join(args.outdir,'best_result.txt.bk'+tmstamp))
                        if exists(join(args.outdir,'plots')):
                            rename(join(args.outdir,'plots'),join(joint(args.outdir,'plots_bk_')+tmstamp))
                print("Creating output file using current max_round beam")
                beamf=join(args.outdir,'beam_{}.p'.format(args.max_round-1))
                curr_beam=OrderedDict(cPickle.load(open(beamf, "rb" )))
                minf=join(args.outdir,'beam_{}'.format(args.max_round-1))
                with open(minf,'r') as f:
                    curr_min=int(f.readline().split('\t')[1])
            else:
                if exists(join(args.outdir,'best_result.txt')):
                    score=float(np.loadtxt(join(args.outdir,'best_result.txt'),dtype='S')[1])
                    cut=len(np.loadtxt(join(args.outdir,'best_result.txt'),dtype='S')[0].split('_'))
                    if score>=args.coverage_cutoff and cut>=args.lower_bound:
                        print("Optimization already reached")
                        sys.exit(0)
                    tmstamp=datetime.now().strftime("%Y%m%d%H%M%S")
                    rename(join(args.outdir,'best_result.txt'),join(args.outdir,'best_result.txt.bk'+tmstamp))
                    if exists(join(args.outdir,'plots')):
                        rename(join(args.outdir,'plots'),join(join(args.outdir,'plots_bk_')+tmstamp))
                beamf=join(args.outdir,'beam_{}.p'.format(max(bnum)))
                curr_beam=OrderedDict(cPickle.load(open(beamf, "rb" )))
                minf=join(args.outdir,'beam_{}'.format(max(bnum)))
                with open(minf,'r') as f:
                    curr_min=int(f.readline().split('\t')[1])
        else:
            curr_beam={}
            curr_min=0

    objectives=pd.read_pickle('AllEpitopeFeatures.pkl')
    selfp=np.loadtxt('self_pept.txt',dtype='S')

    def process_predict(pred_file,corre_file,args,objectives,selfp,basement=None,basefile=None,mtype='mhc1'):
        print "loading file",pred_file
        prediction=pd.read_pickle(pred_file)
        if len(prediction.columns.values[0])==2:
            prediction=prediction.droplevel('loci', axis=1)
        if not args.skippre:
            if args.binary_cutoff:
                print("Binarizing predictions with threshold %f" % args.binary_cutoff)
                prediction=(prediction>=args.binary_cutoff).applymap(lambda x:int(x))
        
        if corre_file:
            print("loading correction files from %s" % corre_file)
            correction=pd.read_pickle(corre_file)
            correction=prediction.loc[correction.index,correction.columns]*correction
        if not args.skippre:
            prediction[objectives.loc[prediction.index]['glyco_probs']>args.glyco_cutoff]=0.0
            prediction[objectives.loc[prediction.index]['crosses_cleavage']>0]=0.0
        if corre_file:
            prediction.loc[correction.index,correction.columns]=correction
        prediction.drop(selfp,errors='ignore',inplace=True)
        if basement:
            if basefile:
                print "loading base design from", basefile
                bpeptides=np.loadtxt(basefile,dtype='S')
                basetable=prediction.loc[bpeptides]
            else:
                if '.txt' in basement:
                    print "loading base protein from", basement
                    bpeptides=np.loadtxt(basement,dtype='S')
                    basetable=prediction.loc[bpeptides]
                else:
                    print "calculating base peptides of {}".format(basement)
                    basetable=prediction[objectives.loc[prediction.index]['protein'].apply(lambda x:(x in basement))]
            basecount=basetable.fillna(0.0).sum(axis=0)
            print "basement peptides:",len(basetable)
            basedir=join(args.outdir,'basement_count_{}.pkl'.format(mtype))
            basecount.to_pickle(basedir)
            pc0=PopulationCoverage(basetable,frequency[mtype],regions,args.outdir)
            base_result=pc0.overall_coverage(epitopes=list(basetable.index),lower=args.lower_bound,verbose=True,pre_map=True,return_pre=True,mtype=mtype)
            base_hist=compute_probs(base_result[1])
            base_hist.to_pickle(join(args.outdir,'lower_bound_hist_{}.pkl'.format(mtype)))
            pre_base=[]
            for cntry in ['White','Black','Asians']:
                pre_base.append(base_result[1][cntry][2])
            pre_base=pd.concat(pre_base,axis=0)
            #print pre_base
            basedir=join(args.outdir,'basement_freq_{}.pkl'.format(mtype))
            pre_base.to_pickle(basedir)
            #remove basement peptides
            if '.txt' in basement:
                print "loading base sequences to remove from", basement
                rmpeptides=np.loadtxt(basement,dtype='S')
                prediction=prediction.drop(rmpeptides,errors='ignore')
            else:
                prediction=prediction[objectives.loc[prediction.index]['protein'].apply(lambda x:(x not in basement))]
        else:
            basecount=None
            base_hist=None
            pre_base=None
        if args.protein:
            prediction=prediction[objectives.loc[prediction.index]['protein'].apply(lambda x:(x in args.protein))]
        prediction=prediction[objectives.loc[prediction.index]['perc_mutated']<=args.mutation_cutoff]
        prediction=prediction[prediction.sum(axis=1)>0]
        print("loaded %d peptides" % len(prediction))
        return prediction,basecount,base_hist,pre_base

    prediction={}
    basecount={}
    base_hist={}
    pre_base={}

    if args.correction:
        corr={'mhc1':'Correction-mhc1.pkl','mhc2':'Correction-mhc2.pkl'}
    else:
        corr={'mhc1':None,'mhc2':None}

    if args.basement or args.basefile:
        for mhc in ['mhc1','mhc2']:
            prediction[mhc],basecount[mhc],base_hist[mhc],pre_base[mhc]=process_predict(pred_files[mhc],corr[mhc],args,objectives,selfp,\
                                                basement=args.basement.format(mhc),basefile=args.basefile.format(mhc),mtype=mhc)
        basedir=join(args.outdir,'basement_freq_{}.pkl') 
        base_hist=(base_hist['mhc1']+base_hist['mhc2'])/2.0   
        print base_hist
    else:
        for mhc in ['mhc1','mhc2']:
            prediction[mhc],basecount[mhc],base_hist[mhc],pre_base[mhc]=process_predict(pred_files[mhc],corr[mhc],args,objectives,selfp,mtype=mhc)
        basedir=None
        base_hist=None

    if not exists(join(args.outdir,'plots')):
        makedirs(join(args.outdir,'plots'))
    candidate_list=pd.read_pickle('preprocess_len25_all_adapt.pkl')
    candidate_list=candidate_list[candidate_list['start_pos']%args.diversity==0]
    candidate_list=candidate_list[(candidate_list['coverage_mhc2']+candidate_list['coverage_mhc1'])>0]
    if args.basement:
        if '.txt' in args.basement:
            print "loading base sequences to remove from", args.basement.format('mhc2')
            rmpeptides=np.loadtxt(args.basement.format('mhc2'),dtype='S')
            candidate_list=candidate_list.drop(rmpeptides,errors='ignore')
        else:
            candidate_list=candidate_list[candidate_list['protein'].apply(lambda x:(x not in args.basement))]

    print 'selecting from length25 peptides:',len(candidate_list)

    #(self,input_epitope,hap_map,country_list,outdir,candidates,pre_map=None,base_counting=None):
    pc=PopulationCoverage2(prediction,frequency,regions,args.outdir,candidate_list,base_counting=basedir)
    if not args.restart:
        #prediction.to_pickle(join(args.outdir,'processed_prediction.pkl'))
        print "calculating maximum coverage with lower bound: %d" % args.lower_bound
        #epitopes_long,lower=0,verbose=False,pre_map=None,return_pre=False):
        max_result=pc.overall_multi(list(candidate_list.index),lower=args.lower_bound,verbose=True,pre_map=True)
        #save_detail(join(args.outdir,'upper_bound.txt'),max_result,has_overall=True)
        #plot_hist(max_result[1],' all peptide',join(args.outdir,'plots','maximum_hist.png'))
        max_hist1=compute_probs(max_result[1]['mhc1'][1])
        max_hist2=compute_probs(max_result[1]['mhc2'][1])
        max_hist=(max_hist1+max_hist2)/2.0
        max_hist.to_pickle(join(args.outdir,'upper_bound_hist.pkl'))
        max_hist1.to_pickle(join(args.outdir,'upper_bound_hist_mhc1.pkl'))
        max_hist2.to_pickle(join(args.outdir,'upper_bound_hist_mhc2.pkl'))
    else:
        max_hist=pd.read_pickle(join(args.outdir,'upper_bound_hist.pkl'))
    print max_hist

    best_solution,detail,beam_history=pc.beam_search_parallel(beam_size=args.beam_size,cutoff=args.coverage_cutoff,min_cutoff=args.lower_bound,\
                    max_round=args.max_round,curr_beam=curr_beam,curr_min=curr_min,nworker=args.nworker, bs=args.batchsize,diverse_cut=args.diversity,augment=(not args.basement is None),upper=max_hist,lower=base_hist)
    with open(join(args.outdir,'best_result.txt'),'w') as f:
        f.writelines('%s\t%f' % best_solution[0],best_solution[1][0])
    #save_detail(join(args.outdir,'best_detail.txt'),detail)

    # print('plotting for peptide set:')
    # print(best_solution[0].split('_'))
    # if best_solution[1]<args.coverage_cutoff:
    #     np.savetxt(join(args.outdir,'failed'),[best_solution[1],args.coverage_cutoff],fmt='%s')
