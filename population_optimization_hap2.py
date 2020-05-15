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
        return OrderedDict([(k,v) for v,k in sorted(self.heap,reverse=True)])

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
    pre=pd.read_pickle('preprocess_haplotype_new{}.pkl'.format(args.type.split('_')[0][-1]))
    premap={}
    for c in ['White','Black','Asians']:
        premap[c]=pre[pre['country']==c].drop(labels='country',axis=1)
    
class PopulationCoverage(object):
    def __init__(self,input_epitope,hap_map,country_list,outdir,pre_map=None):
        self.input_epitope_affinity=input_epitope
        self.hap_map=hap_map
        self.outdir=outdir
        if pre_map:
            #self.pre_map=pre_map
            self.set_country(country_list,precompute=False)
        else:
            self.set_country(country_list,precompute=False)

    def beam_search_parallel2(self,beam_size=20,cutoff=0.9,min_cutoff=5,\
                            max_round=20,curr_beam={},curr_min=0,nworker=20,bs=50,diverse_cut=3):
        print 'Using %d workers' % nworker
        current_max_coverage=curr_beam.items()[0][1] if len(curr_beam)>0 else 0.0
        beams=[]
        curr_length=len(curr_beam.items()[0][0].split('_')) if len(curr_beam)>0 else 0
        pool=pp._ProcessPool(processes=nworker,initializer=initializer,initargs=())
        outdir=join(self.outdir,'plots')
        print 'current beam: ',curr_beam.items()
        print 'current lower bound: ',curr_min
        while (current_max_coverage<cutoff or curr_min<min_cutoff) and curr_length<max_round:
            print 'beamsearch round ',curr_length
            t0 = time.time()
            self.next_beam=KthLargest(k=beam_size)
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x.split('_') for x in curr_beam.keys()]
            all_args=[]
            print 'current candidates:',curr_candidates
            test_epitopes=[x+[y] for x in curr_candidates for y in self.input_epitope_affinity.index if diff1d(x,y,cut=diverse_cut)]
            for i in range(0,len(test_epitopes),bs):
                all_args.append((copy(test_epitopes[i:min(i+bs,len(test_epitopes))]),curr_min))
            print len(test_epitopes),len(all_args)
            del test_epitopes
            r=pool.map_async(self.batch_scan,all_args,callback=self.batch_update)
            r.get()
            curr_beam=self.next_beam.return_dict()
            current_max_coverage=curr_beam.items()[0][1]
            beams.append(curr_beam)
            # check histogram here
            curr_length=len(curr_beam.items()[0][0].split('_'))
            # for k,e in enumerate(curr_beam.items()):
            #     cover, hist_detail= self.overall_coverage(epitopes=e[0].split('_'), lower=curr_min,pre_map=True,verbose=True)
            #     plot_hist(hist_detail,' {:.2f} for lb={} #pept={}'.format(cover,curr_min,curr_length),join(outdir,'{}_{}_hist_{}.png'.format(curr_length,curr_min,k)))
            current_median=np.median([v[1] for v in curr_beam.items()])
            old_min=curr_min
            if current_median > cutoff:
                print 'median min_coverage of {} reached {}, raising min_coverage'.format(curr_min,current_median)
                curr_min+=1
            print 'current beam: ',curr_beam.items()
            print 'current lower bound: ',curr_min
            save_pickle(join(self.outdir,'beam_'+str(curr_length-1)+'.p'),curr_beam.items())
            save_beams(join(self.outdir,'beam_'+str(curr_length-1)),curr_beam,curr_min,old_min)
            t1 = time.time()
            print('time passed: {}'.format(t1-t0))

        pool.close()
        print 'Coverage cutoff reached, final solution:',curr_beam.items()[0]
        print 'Per region details:'
        details=self.overall_coverage(epitopes=curr_beam.items()[0][0].split('_'),lower=min_cutoff,pre_map=True,verbose=True)
        plot_hist(details[1],'final {:.2f} for lb={} #pept={}'.format(details[0],min_cutoff,curr_length),join(outdir,'final_hist.png'))
        return curr_beam.items()[0],details,beams

    def beam_search_parallel(self,beam_size=20,cutoff=0.9,min_cutoff=5,\
                            max_round=20,curr_beam={},curr_min=0,nworker=20,bs=50,diverse_cut=3):
        print 'Using %d workers' % nworker
        current_max_coverage=curr_beam.items()[0][1] if len(curr_beam)>0 else 0.0
        beams=[]
        curr_length=len(curr_beam.items()[0][0].split('_')) if len(curr_beam)>0 else 0
        outdir=join(self.outdir,'plots')
        print 'current beam: ',curr_beam.items()
        print 'current lower bound: ',curr_min
        if curr_min>0:
                iter_cutoff=cutoff
        else:
            if args.initial_cut:
                iter_cutoff=max(args.initial_cut,cutoff)
            else:
                if args.type.split('_')[0][-1]=='1':
                    iter_cutoff=max(0.97,cutoff)
                else:
                    iter_cutoff=max(0.93,cutoff)
        while (current_max_coverage<iter_cutoff or curr_min<min_cutoff) and curr_length<max_round:
            pool=pp._ProcessPool(processes=nworker,initializer=initializer,initargs=())
            print 'beamsearch round ',curr_length, 'cutoff', iter_cutoff
            t0 = time.time()
            self.next_beam=KthLargest(k=beam_size)
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x.split('_') for x in curr_beam.keys()]
            all_args=[]
            print 'current candidates:',curr_candidates
            test_epitopes=[]
            seen_curr={}
            for c_keys in curr_candidates:
#                 counting=self.input_epitope_affinity.loc[c_keys].fillna(0.0).sum(axis=0)
#                 curr_key=','.join(str(x) for x in counting.values)
                ids=[x for x in self.input_epitope_affinity.index if diff1d(c_keys,x,cut=diverse_cut)]
                part=self.input_epitope_affinity.loc[ids].copy()
                part['code']=part.apply(lambda x:''.join([str(k) for k in x]),axis=1)
                for name, group in part.groupby('code'):
                    if len(group)>2 and curr_length>0:
                        group=group.sample(n=2)
                    test_epitopes.append((c_keys,group.index))
            for divide in range(10):
                bs=len(test_epitopes)//((divide+1)*nworker)+int(len(test_epitopes)%((divide+1)*nworker)>0)
                if bs<60:
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
            current_max_coverage=curr_beam.items()[0][1]
            beams.append(curr_beam)
            # check histogram here
            curr_length=len(curr_beam.items()[0][0].split('_'))
            # for k,e in enumerate(curr_beam.items()):
            #     cover, hist_detail= self.overall_coverage(epitopes=e[0].split('_'), lower=curr_min,pre_map=True,verbose=True)
            #     plot_hist(hist_detail,' {:.2f} for lb={} #pept={}'.format(cover,curr_min,curr_length),join(outdir,'{}_{}_hist_{}.png'.format(curr_length,curr_min,k)))
            current_median=np.median([v[1] for v in curr_beam.items()])
            old_min=curr_min
            if current_median > iter_cutoff:
                print 'median min_coverage of {} reached {}, raising min_coverage'.format(curr_min,current_median)
                curr_min+=1
            print 'current beam: ',curr_beam.items()
            print 'current lower bound: ',curr_min
            save_pickle(join(self.outdir,'beam_'+str(curr_length-1)+'.p'),curr_beam.items())
            save_beams(join(self.outdir,'beam_'+str(curr_length-1)),curr_beam,curr_min,old_min)
            t1 = time.time()
            print('time passed: {}'.format(t1-t0))
            pool.close()
            if curr_min>0:
                iter_cutoff=cutoff

        print 'Coverage cutoff reached, final solution:',curr_beam.items()[0]
        print 'Per region details:'
        details=self.overall_coverage(epitopes=curr_beam.items()[0][0].split('_'),lower=min_cutoff,pre_map=True,verbose=True)
        plot_hist(details[1],'final {:.2f} for lb={} #pept={}'.format(details[0],min_cutoff,curr_length),join(outdir,'final_hist.png'))
        return curr_beam.items()[0],details,beams

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
        for test_setb in test_sets:
            old,new=test_setb
            for i,new_seq in enumerate(new):
                test_set=old+[new_seq]
                test_set=sorted(test_set)
                key='_'.join(test_set)
                if i==0:
                    coverage=self.overall_coverage(test_set,lower_bound,verbose=False)
                result.append((coverage,key))
        for item in result:
                local_beam.add(item)
        return [(x[1],x[0]) for x in local_beam.return_dict().items()]

    def batch_update(self,results):
        print 'called update'
        print len(results[0]),len(results)
        for rs in results:
            for item in rs:
                self.next_beam.add(item)

    def overall_coverage(self,epitopes,lower=0,verbose=True,pre_map=None):
        country_coverage=[]
        country_detail={}
        if isinstance(epitopes,str):
            epitopes=[epitopes]
        counting=self.input_epitope_affinity.loc[epitopes].fillna(0.0).sum(axis=0)
        for country in self.country_list:
            if pre_map:
                result=self._compute_coverage(country,counting,lower,verbose=verbose)
            else:
                result=self.compute_coverage(country,counting,lower,verbose=verbose)
            country_coverage.append(result[0])
            if verbose:
                country_detail[country]=result
        if verbose:
            return np.mean(country_coverage),country_detail
        else:
            return np.mean(country_coverage)
    
    def compute_coverage(self,country,counting,lower,verbose=True,truncate=True):
        global premap
        new_pre3=premap[country].copy()
        #new_pre3['count']=0
        valid=set(counting[counting>0].index)&self.country_allele[country]
        prob,hist=self._compute_hist(valid,new_pre3,counting,lower,args.type)
        if verbose:
            return prob,hist
        else:
            return [prob]

    def _compute_coverage(self,country,counting,lower,verbose=True,truncate=True):
        pre=pd.read_pickle('preprocess_haplotype_new{}.pkl'.format(args.type.split('_')[0][-1]))
        new_pre3=pre[pre['country']==country]#.drop(labels='country',axis=1)
        #new_pre3['count']=0
        valid=set(counting[counting>0].index)&self.country_allele[country]
        #print 'searching on %d alleles' % len(valid)
        prob,hist=self._compute_hist(valid,new_pre3,counting,lower,args.type)
        if verbose:
            return prob,hist
        else:
            return [prob]

    def _compute_hist(self,valid,new_pre3,counting,lower,type):
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
    parser.add_argument("-ic","--initcutoff",dest='initial_cut', type=float,
                        help="Target coverage lower bond when stop beam search")
    parser.add_argument("-lo","--lower",dest='lower_bound', type=int, default=3,
                        help="Target coverage lower bond #peptide when stop beam search")
    parser.add_argument("-r","--maxround",dest='max_round', type=int, default=35,
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
    parser.add_argument("-d","--diversity", type=int, default=3,
                        help="Distance for removing similar sliding windows")
    parser.add_argument("--unroll", action='store_false',
                        help="path for results.")
    parser.add_argument("--downsamp", action='store_true',
                        help="path for results.")
    parser.add_argument("--restart", action='store_true',
                        help="If restarting from previous results or not")
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
            f.writelines('%s\t%f\n' % (seq,val))

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
        else:
            frequency=pd.read_pickle('haplotype_frequency_marry2.pkl')
    else:
        print("prediction output type unknown")
        sys.exit(0)

    if not exists(args.regions):
        print("Region file not exist: %s" % args.regions)
        sys.exit(0)
    regions=['White','Black','Asians']#np.loadtxt(args.regions,dtype='S',delimiter='\n')

    if args.method not in ['puffin','deepligand','mhcflurry','netmhc','test','average','mean','netmhcii-3.2','netmhcii-4.0']:
        print("prediction method unknown")
        sys.exit(0)
    if args.prediction not in ['pred_prob','ic50','pred_affinity','rank']:
        print("prediction output type unknown")
        sys.exit(0)

    pred_file='_'.join([args.type, args.method, args.prediction, 'pivot.pkl.gz'])
    if not exists(pred_file):
        print("prediction file not exist %s" % pred_file)
        sys.exit(0)
    print('reading prediction file: %s' % pred_file)
    prediction=pd.read_pickle(pred_file).droplevel('loci', axis=1)
    prediction_original=prediction.copy()
    if args.binary_cutoff:
        print("Binarizing predictions with threshold %f" % args.binary_cutoff)
        # if args.binary_cutoff>1:
        #     prediction=(prediction<args.binary_cutoff).applymap(lambda x:int(x))
        # else:
        prediction=(prediction>=args.binary_cutoff).applymap(lambda x:int(x))
        # idx = pd.IndexSlice
        # prediction.loc[:,idx[:,'unknown']]=0.0
    elif args.truncate_cutoff:
        print("Truncating predictions with threshold %f" % args.truncate_cutoff)
        prediction=prediction.applymap(lambda x:0 if x<args.truncate_cutoff else x)
    
    if prediction.max().max()>1 or prediction.min().min()<0:
        print("prediction values exceeded normal range (%f,%f)" % (prediction.min().min(),prediction.max().max()))
        sys.exit(0)
    
    # if len(prediction.columns)>len(frequency.columns):
    #     prediction=prediction[frequency.columns.values]
    #     prediction.columns=frequency.columns
    # elif len(prediction.columns)<len(frequency.columns):
    #     print("mismatch between allele numbers (%d,%d)" % (len(prediction.columns),len(frequency.columns)))
    #     sys.exit(0)
    # print("loaded %d alleles" % (len(frequency.columns)-3))

    objectives=pd.read_pickle('AllEpitopeFeatures.pkl')
    selfp=np.loadtxt('self_pept.txt',dtype='S')
    prediction.drop(selfp,errors='ignore',inplace=True)
    prediction=prediction[objectives.loc[prediction.index]['glyco_probs']<=args.glyco_cutoff]
    if args.protein:
        prediction=prediction[objectives.loc[prediction.index]['protein'].apply(lambda x:(x in args.protein))]
    prediction=prediction[objectives.loc[prediction.index]['crosses_cleavage']==0]
    prediction=prediction[objectives.loc[prediction.index]['perc_mutated']<=args.mutation_cutoff]
    prediction=prediction[prediction.sum(axis=1)>0]
    if args.downsamp:
        prediction=prediction.sample(n=200)
    print("loaded %d peptides" % len(prediction))
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
                    if score>=args.coverage_cutoff:
                        print("Optimization already reached")
                        sys.exit(0)
                    tmstamp=datetime.now().strftime("%Y%m%d%H%M%S")
                    rename(join(args.outdir,'best_result.txt'),join(args.outdir,'best_result.txt.bk'+tmstamp))
                    if exists(join(args.outdir,'plots')):
                        rename(join(args.outdir,'plots'),join(joint(args.outdir,'plots_bk_')+tmstamp))
                beamf=join(args.outdir,'beam_{}.p'.format(max(bnum)))
                curr_beam=OrderedDict(cPickle.load(open(beamf, "rb" )))
                minf=join(args.outdir,'beam_{}'.format(max(bnum)))
                with open(minf,'r') as f:
                    curr_min=int(f.readline().split('\t')[1])
        else:
            curr_beam={}
            curr_min=0

    if not exists(join(args.outdir,'plots')):
        makedirs(join(args.outdir,'plots'))

    pc=PopulationCoverage(prediction,frequency,regions,args.outdir)
    if not args.restart:
        prediction.to_pickle(join(args.outdir,'processed_prediction.pkl'))
        print "calculating maximum coverage with lower bound: %d" % args.lower_bound
        max_result=pc.overall_coverage(epitopes=list(prediction.index),lower=args.lower_bound,verbose=True,pre_map=True)
        save_detail(join(args.outdir,'upper_bound.txt'),max_result,has_overall=True)
        plot_hist(max_result[1],' all peptide',join(args.outdir,'plots','maximum_hist.png'))
    if args.nworker:
        best_solution,detail,beam_history=pc.beam_search_parallel(beam_size=args.beam_size,cutoff=args.coverage_cutoff,min_cutoff=args.lower_bound,\
                    max_round=args.max_round,curr_beam=curr_beam,curr_min=curr_min,nworker=args.nworker, bs=args.batchsize,diverse_cut=args.diversity)
    else:
        best_solution,detail,beam_history=pc.beam_search(beam_size=args.beam_size,cutoff=args.coverage_cutoff,max_round=args.max_round,roll=args.unroll,curr_beam=curr_beam,curr_min=curr_min)
    with open(join(args.outdir,'best_result.txt'),'w') as f:
        f.writelines('%s\t%f' % best_solution)
    save_detail(join(args.outdir,'best_detail.txt'),detail)

    print('plotting for peptide set:')
    print(best_solution[0].split('_'))
    #plot_detail(prediction,best_solution[0].split('_'),pc,join(args.outdir,'plots'))
    #plot_detail(prediction_original,best_solution[0].split('_'),pc,join(args.outdir,'plots'),'_raw',norm=(args.prediction!='ic50'))
    plot_detail3(prediction, best_solution[0].split('_'), pc,join(args.outdir,'plots'),lgd=len(best_solution[0].split('_'))<20)
    plot_detail3(prediction_original, best_solution[0].split('_'), pc, join(args.outdir,'plots'),suffix='_raw',lgd=len(best_solution[0].split('_'))<20)
    if best_solution[1]<args.coverage_cutoff:
        np.savetxt(join(args.outdir,'failed'),[best_solution[1],args.coverage_cutoff],fmt='%s')
