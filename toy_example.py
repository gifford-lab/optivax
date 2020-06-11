import pandas as pd
import numpy as np
import seaborn as sns
from collections import OrderedDict
import heapq
from tqdm import tqdm
import pylab,random,cPickle
from copy import copy
import time
from os.path import dirname, basename,join,exists
from os import makedirs,system,listdir,rmdir
import multiprocessing as mp
#import cPickle
import pathos.pools as pp
from argparse import ArgumentParser
#import faulthandler; faulthandler.enable()

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
    def __init__(self,input_epitope,hap_map):
        self.input_epitope_affinity=input_epitope
        self.hap_map=hap_map

    def single_search(self,beam_size=1,cutoff=0.9,min_cutoff=3,max_round=20,nworker=20,bs=100,look=True):
        current_max_coverage=0.0
        beams=[]
        curr_beam=[]
        cnt=0
        pool=pp._ProcessPool(processes=nworker,initargs=())
        curr_min=1
        up=float('inf')
        while cnt<max_round:
            t0 = time.time()
            #print 'round ',cnt
            self.next_beam=[]
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x[1].split('_') for x in curr_beam]#.keys()]
            all_args=[]
            test_epitopes=[x+[y] for x in curr_candidates for y in self.input_epitope_affinity.index if not y in x]
            for i in range(0,len(test_epitopes),bs):
                all_args.append((copy(test_epitopes[i:min(i+bs,len(test_epitopes))]),curr_min))
            del test_epitopes
            self.progress=0
            r=pool.map_async(self.batch_scan,all_args,callback=self.batch_update)
            r.get()
            scores=[x[0][x[0]['count']>=curr_min]['freq'].sum() for x in self.next_beam]
            topidx=np.argwhere(scores == np.amax(scores)).flatten()
            if look:
                tb_min=curr_min
                while len(topidx)>1 and tb_min<=min_cutoff:
                    self.next_beam=self.next_beam[list(topidx)]
                    tb_min+=1
                    scores=[x[0][x[0]['count']>=tb_min]['freq'].sum() for x in self.next_beam]
                    topidx=np.argwhere(scores == np.amax(scores)).flatten()
            curr_beam=[self.next_beam[np.random.choice(topidx,1)[0]]]
            beams.append(curr_beam[0])
            h=curr_beam[0][0]
            current_max_coverage=h[h['count']>=curr_min]['freq'].sum()
            if current_max_coverage >= cutoff:
                if curr_min==min_cutoff:
                    break
                else:
                    curr_min+=1
            cnt+=1
            t1 = time.time()
        pool.close()
        return curr_beam,beams
    
    def greedy_search(self,beam_size=1,cutoff=0.9,min_cutoff=3,max_round=20,nworker=20,bs=100):
        current_max_coverage=0.0
        beams=[]
        curr_beam=[]
        cnt=0
        pool=pp._ProcessPool(processes=nworker,initializer=initializer,initargs=())
        curr_min=min_cutoff
        up=float('inf')
        while current_max_coverage<cutoff and cnt<max_round:
            t0 = time.time()
            self.next_beam=[]
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x[1].split('_') for x in curr_beam]#.keys()]
            all_args=[]
            test_epitopes=[x+[y] for x in curr_candidates for y in self.input_epitope_affinity.index if not y in x]
            for i in range(0,len(test_epitopes),bs):
                all_args.append((copy(test_epitopes[i:min(i+bs,len(test_epitopes))]),curr_min))
            del test_epitopes
            r=pool.map_async(self.batch_scan,all_args,callback=self.batch_update)
            r.get()
            scores=[x[0][x[0]['count']>=curr_min]['freq'].sum() for x in self.next_beam]
            topidx=np.argwhere(scores == np.amax(scores)).flatten()
            curr_beam=[self.next_beam[np.random.choice(topidx,1)[0]]]
            h=curr_beam[0][0]
            current_max_coverage=h[h['count']>=curr_min]['freq'].sum()
            beams.append((current_max_coverage,curr_beam[0][1]))
            cnt+=1
            t1 = time.time()
        pool.close()
        return curr_beam,beams
    
    def exhaust_search(self,beam_size=1,cutoff=0.9,min_cutoff=3,max_round=20,nworker=20,bs=100):
        current_max_coverage=0.0
        beams=[]
        curr_beam=[]
        cnt=0
        pool=pp._ProcessPool(processes=nworker,initializer=initializer,initargs=())
        curr_min=min_cutoff
        up=float('inf')
        while current_max_coverage<cutoff and cnt<max_round:
            t0 = time.time()
            print('round ',cnt)
            self.next_beam=[]
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x[1].split('_') for x in curr_beam]#.keys()]
            all_args=[]
            test_epitopes=[x+[y] for x in curr_candidates for y in self.input_epitope_affinity.index if not y in x]
            for i in range(0,len(test_epitopes),bs):
                all_args.append((copy(test_epitopes[i:min(i+bs,len(test_epitopes))]),curr_min))
            print(len(test_epitopes),len(all_args))
            del test_epitopes
            r=pool.map_async(self.batch_scan,all_args,callback=self.batch_update)
            r.get()
            scores=[x[0] for x in self.next_beam]
            curr_beam=self.next_beam
            current_max_coverage=np.max(scores)
            beams.append((current_max_coverage,self.next_beam[np.argmax(scores)]))
            cnt+=1
            t1 = time.time()
            print('time passed: {}'.format(t1-t0))
        pool.close()
        return curr_beam,beams

    def batch_scan(self,inputs):
        test_sets,lower_bound = inputs
        result=[]
        for test_set in test_sets:
            test_set=sorted(test_set)
            key='_'.join(test_set)
            hist=self.overall_coverage(test_set,verbose=False)
            p=hist[hist['count']>=lower_bound]['freq'].sum()
            result.append((p,key))
        return result

    def batch_update(self,results):
        print('called update')
        for rs in results:
            for item in rs:
                self.next_beam.append(item)
        #self.next_beam=np.asarray(self.next_beam)

    def overall_coverage(self,epitopes,verbose=True):
        if isinstance(epitopes,str):
            epitopes=[epitopes]
        counting=self.input_epitope_affinity.loc[epitopes].fillna(0.0).sum(axis=0)
        weight=self.hap_map
        new_pre3=pd.concat([counting.to_frame('count'),weight.to_frame('freq')],axis=1)
        hist=new_pre3.groupby('count').sum().reset_index()
        return hist

class PopulationCoverage4(object):
    def __init__(self,input_epitope,hap_map):
        self.input_epitope_affinity=input_epitope
        self.hap_map=hap_map

    def single_search(self,beam_size=1,cutoff=0.9,min_cutoff=3,max_round=20,nworker=20,bs=100,look=True):
        current_max_coverage=0.0
        beams=[]
        curr_beam=[]
        cnt=0
        pool=pp._ProcessPool(processes=nworker,initializer=initializer,initargs=())
        curr_min=1
        up=float('inf')
        while cnt<max_round:
            t0 = time.time()
            #print 'round ',cnt
            self.next_beam=[]
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x[1].split('_') for x in curr_beam]#.keys()]
            all_args=[]
            test_epitopes=[x+[y] for x in curr_candidates for y in self.input_epitope_affinity.index if not y in x]
            for i in range(0,len(test_epitopes),bs):
                all_args.append((copy(test_epitopes[i:min(i+bs,len(test_epitopes))]),curr_min))
            del test_epitopes
            self.progress=0
            r=pool.map_async(self.batch_scan,all_args,callback=self.batch_update)
            r.get()
            scores=[x[0][x[0]['count']>=curr_min]['freq'].sum() for x in self.next_beam]
            topidx=np.argwhere(scores == np.amax(scores)).flatten()
            if look:
                tb_min=curr_min
                while len(topidx)>1 and tb_min<=min_cutoff:
                    self.next_beam=self.next_beam[list(topidx)]
                    tb_min+=1
                    scores=[x[0][x[0]['count']>=tb_min]['freq'].sum() for x in self.next_beam]
                    topidx=np.argwhere(scores == np.amax(scores)).flatten()
            curr_beam=[self.next_beam[np.random.choice(topidx,1)[0]]]
            beams.append(curr_beam[0])
            h=curr_beam[0][0]
            current_max_coverage=h[h['count']>=curr_min]['freq'].sum()
            if current_max_coverage >= cutoff:
                    curr_min+=1
            cnt+=1
            t1 = time.time()
        pool.close()
        return curr_beam,beams
    
    def greedy_search(self,beam_size=1,cutoff=0.9,min_cutoff=3,max_round=20,nworker=20,bs=100):
        current_max_coverage=0.0
        beams=[]
        curr_beam=[]
        cnt=0
        pool=pp._ProcessPool(processes=nworker,initializer=initializer,initargs=())
        curr_min=min_cutoff
        up=float('inf')
        while cnt<max_round:
            t0 = time.time()
            #print 'round ',cnt
            self.next_beam=[]
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x[1].split('_') for x in curr_beam]#.keys()]
            all_args=[]
            test_epitopes=[x+[y] for x in curr_candidates for y in self.input_epitope_affinity.index if not y in x]
            for i in range(0,len(test_epitopes),bs):
                all_args.append((copy(test_epitopes[i:min(i+bs,len(test_epitopes))]),curr_min))
            del test_epitopes
            r=pool.map_async(self.batch_scan,all_args,callback=self.batch_update)
            r.get()
            scores=[x[0][x[0]['count']>=curr_min]['freq'].sum() for x in self.next_beam]
            topidx=np.argwhere(scores == np.amax(scores)).flatten()
            curr_beam=[self.next_beam[np.random.choice(topidx,1)[0]]]
            h=curr_beam[0][0]
            current_max_coverage=h[h['count']>=curr_min]['freq'].sum()
            beams.append((current_max_coverage,curr_beam[0][1]))
            cnt+=1
            t1 = time.time()
        pool.close()
        return curr_beam,beams
    
    def exhaust_search(self,beam_size=1,cutoff=0.9,min_cutoff=3,max_round=20,nworker=20,bs=100):
        current_max_coverage=0.0
        beams=[]
        curr_beam=[]
        cnt=0
        pool=pp._ProcessPool(processes=nworker,initializer=initializer,initargs=())
        curr_min=min_cutoff
        up=float('inf')
        while cnt<max_round:
            t0 = time.time()
            print('round ',cnt)
            self.next_beam=[]
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x[1].split('_') for x in curr_beam]#.keys()]
            all_args=[]
            test_epitopes=[x+[y] for x in curr_candidates for y in self.input_epitope_affinity.index if not y in x]
            for i in range(0,len(test_epitopes),bs):
                all_args.append((copy(test_epitopes[i:min(i+bs,len(test_epitopes))]),curr_min))
            print(len(test_epitopes),len(all_args))
            del test_epitopes
            r=pool.map_async(self.batch_scan,all_args,callback=self.batch_update)
            r.get()
            scores=[x[0] for x in self.next_beam]
            curr_beam=self.next_beam
            current_max_coverage=np.max(scores)
            beams.append(self.next_beam[np.argmax(scores)])
            cnt+=1
            t1 = time.time()
            print('time passed: {}'.format(t1-t0))
        pool.close()
        print('Coverage cutoff reached, final solution:',beams[-1])
        return curr_beam,beams

    def batch_scan(self,inputs):
        #print 'called batch'
        test_sets,lower_bound = inputs
        result=[]
        for test_set in test_sets:
            test_set=sorted(test_set)
            key='_'.join(test_set)
            hist=self.overall_coverage(test_set,verbose=False)
            p=hist[hist['count']>=lower_bound]['freq'].sum()
            result.append((p,key))
        return result

    def batch_update(self,results):
        print('called update')
        for rs in results:
            for item in rs:
                self.next_beam.append(item)#.add(item)
        #self.next_beam=np.asarray(self.next_beam)

    def overall_coverage(self,epitopes,verbose=True):
        if isinstance(epitopes,str):
            epitopes=[epitopes]
        counting=self.input_epitope_affinity.loc[epitopes].fillna(0.0).sum(axis=0)
        weight=self.hap_map
        new_pre3=pd.concat([counting.to_frame('count'),weight.to_frame('freq')],axis=1)
        hist=new_pre3.groupby('count').sum().reset_index()
        return hist

class PopulationCoverage5(object):
    def __init__(self,input_epitope,hap_map):
        self.input_epitope_affinity=input_epitope
        self.hap_map=hap_map    
    def exhaust_search(self,beam_size=1,cutoff=0.9,min_cutoff=3,max_round=20,nworker=20,bs=100):
        current_max_coverage=0.0
        beams=[]
        curr_beam=[]
        cnt=0
        pool=pp._ProcessPool(processes=nworker,initializer=initializer,initargs=())
        curr_min=min_cutoff
        up=float('inf')
        while cnt<max_round:
            t0 = time.time()
            print 'round ',cnt
            self.next_beam=[]
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x[1].split('_') for x in curr_beam]#.keys()]
            all_args=[]
            test_epitopes=[x+[y] for x in curr_candidates for y in self.input_epitope_affinity.index if not y in x]
            for i in range(0,len(test_epitopes),bs):
                all_args.append((copy(test_epitopes[i:min(i+bs,len(test_epitopes))]),curr_min))
            print len(test_epitopes),len(all_args)
            del test_epitopes
            r=pool.map_async(self.batch_scan,all_args,callback=self.batch_update)
            r.get()
            scores=[x[0] for x in self.next_beam]#[x[0]['count']>=curr_min]['freq'].sum()
            #print scores
            curr_beam=self.next_beam
            current_max_coverage=np.max(scores)
            beams.append(self.next_beam[np.argmax(scores)])
            cnt+=1
            t1 = time.time()
            print('time passed: {}'.format(t1-t0))
        pool.close()
        return curr_beam,beams
    
    def exhaust_search_cover(self,beam_size=1,cutoff=0.9,min_cutoff=3,max_round=20,nworker=20,bs=100):
        current_max_coverage=0.0
        beams=[]
        curr_beam=[]
        cnt=0
        pool=pp._ProcessPool(processes=nworker,initializer=initializer,initargs=())
        curr_min=min_cutoff
        up=float('inf')
        while current_max_coverage<cutoff and cnt<max_round:
            t0 = time.time()
            print 'round ',cnt
            self.next_beam=[]
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x[1].split('_') for x in curr_beam]#.keys()]
            all_args=[]
            test_epitopes=[x+[y] for x in curr_candidates for y in self.input_epitope_affinity.index if not y in x]
            for i in range(0,len(test_epitopes),bs):
                all_args.append((copy(test_epitopes[i:min(i+bs,len(test_epitopes))]),curr_min))
            #print len(test_epitopes),len(all_args)
            del test_epitopes
            r=pool.map_async(self.batch_scan,all_args,callback=self.batch_update)
            r.get()
            scores=[x[0] for x in self.next_beam]#[x[0][x[0]['count']>=curr_min]['freq'].sum() for x in self.next_beam]
            curr_beam=self.next_beam
            current_max_coverage=np.max(scores)
            beams.append((current_max_coverage,self.next_beam[np.argmax(scores)]))
            cnt+=1
            t1 = time.time()
            print('time passed: {}'.format(t1-t0))
        pool.close()
        return curr_beam,beams

    def batch_scan(self,inputs):
        #print 'called batch'
        test_sets,lower_bound = inputs
        result=[]
        for test_set in test_sets:
            test_set=sorted(test_set)
            key='_'.join(test_set)
            hist=self.overall_coverage(test_set,verbose=False)
            p=hist[hist['count']>=lower_bound]['freq'].sum()
            result.append((p,key))
            #result.append((hist,key))
        return result

    def batch_update(self,results):
        print 'called update'
        for rs in results:
            for item in rs:
                self.next_beam.append(item)#.add(item)
        #self.next_beam=np.asarray(self.next_beam)
        #self.progress+=len(results)
        #print self.progress

    def overall_coverage(self,epitopes,verbose=True):
        if isinstance(epitopes,str):
            epitopes=[epitopes]
        counting=self.input_epitope_affinity.loc[epitopes].fillna(0.0).sum(axis=0)
        weight=self.hap_map
        new_pre3=pd.concat([counting.to_frame('count'),weight.to_frame('freq')],axis=1)
        hist=new_pre3.groupby('count').sum().reset_index()
        return hist

def gen_prob(num_over,num_ele,seed=0,weight=None,cweight=None,upper=3):
    np.random.seed(seed)
    if weight:
        hap=pd.Series(weight)
    else:
        hap=pd.Series([1.0/num_ele]*num_ele)
    hap.index=['x'+str(i) for i in range(num_ele)]
    if cweight:
        cover=np.random.choice(range(upper),size=(num_over,num_ele),p=cweight)
    else:
        cover=np.random.choice(range(upper),size=(num_over,num_ele))
    overlay=pd.DataFrame(cover)
    overlay.columns=hap.index
    overlay.index=['overlay'+str(i) for i in range(num_over)]
    return hap,overlay

def run_compare(seeds,num_over,num_ele,cutoff=5,coverage=0.9,weight=None,cweight=None,max_round=None,upper=3):
    result_pd=[]
    our_pd=[]
    if max_round:
        print "running maximum coverage setup"
        for sd in seeds:
            print sd
            hap,over=gen_prob(num_over,num_ele,seed=sd,upper=upper)
            pc1=PopulationCoverage4(over,hap)
            rs_greed=pc1.greedy_search(cutoff=coverage,min_cutoff=cutoff,max_round=max_round,nworker=20,bs=100)
            for val,key in rs_greed[1]:
                result_pd.append((len(key.split('_')),val,'Greedy'))
            rs_our=pc1.single_search(cutoff=coverage,min_cutoff=cutoff,max_round=max_round,nworker=20,bs=100)
            for hist,key in rs_our[1]:
                for c in range(cutoff):
                    p=hist[hist['count']>=c+1]['freq'].sum()
                    our_pd.append((len(key.split('_')),p,'n={}'.format(c+1)))
                    if c+1==cutoff:
                        result_pd.append((len(key.split('_')),p,'MargiGreedy(Ours)'))
            for i in range(max_round):
                hist=pc1.overall_coverage(over.sample(i+1).index)
                p=hist[hist['count']>=cutoff]['freq'].sum()
                result_pd.append((i+1,p,'Random subset'))
        result_pd=pd.DataFrame(result_pd,columns=['Num. overlays used','n-time coverage','method'])
        our_pd=pd.DataFrame(our_pd,columns=['Num. overlays used','n-time coverage','n'])
        plt.figure(figsize=(4,3))
        ax = sns.lineplot(x='Num. overlays used', y="n-time coverage", hue="method",data=result_pd)
        plt.title("Maximum n-time coverage task (n={})".format(cutoff))
        plt.legend(loc=4,framealpha=0.4)
        plt.figure(figsize=(4,3))
        ax = sns.lineplot(x='Num. overlays used', y="n-time coverage", hue="n",data=our_pd)
        plt.show()
    else:
        print "running set cover setup"
        for sd in seeds:
            print sd
            hap,over=gen_prob(num_over,num_ele,seed=sd,upper=upper)
            pc1=PopulationCoverage3(over,hap)
            for n in range(cutoff):
                print n
                rs_greed=pc1.greedy_search(cutoff=coverage,min_cutoff=n+1,max_round=100,nworker=20,bs=100)
                result_pd.append((n+1,len(rs_greed[0][0][1].split('_')),'Greedy'))
                rs_our=pc1.single_search(cutoff=coverage,min_cutoff=n+1,max_round=100,nworker=20,bs=100)
                result_pd.append((n+1,len(rs_our[0][0][1].split('_')),'MargiGreedy(Ours)'))
        result_pd=pd.DataFrame(result_pd,columns=['n','Num. overlays used','method'])
        plt.figure(figsize=(4,3))
        ax = sns.lineplot(x="n", y="Num. overlays used", hue="method",data=result_pd)
        plt.title("n-time set cover task (0<n<={})".format(cutoff))
        plt.legend(loc=4,framealpha=0.4)
        plt.show()
    return result_pd

def run_compare2(seeds,num_over,num_ele,cutoff=5,coverage=0.9,weight=None,cweight=None,max_round=None,nworker=60,upper=3):
    result_pd=[]
    our_pd=[]
    if max_round:
        print("running maximum coverage setup")
        for sd in seeds:
            print(sd)
            hap,over=gen_prob(num_over,num_ele,seed=sd,upper=upper)
            pc1=PopulationCoverage5(over,hap)
            rs_greed=pc1.exhaust_search(cutoff=coverage,min_cutoff=cutoff,max_round=max_round,nworker=nworker,bs=1000)
            for val,key in rs_greed[1]:
                result_pd.append((len(key.split('_')),val,'True Optimal'))
        result_pd=pd.DataFrame(result_pd,columns=['Num. overlays used','n-time coverage','method'])
    else:
        print("running set cover setup")
        for sd in seeds:
            print(sd)
            hap,over=gen_prob(num_over,num_ele,seed=sd,upper=upper)
            pc1=PopulationCoverage5(over,hap)
            for n in range(cutoff):
                print(n)
                rs_greed=pc1.exhaust_search_cover(cutoff=coverage,min_cutoff=n+1,max_round=100,nworker=nworker,bs=1000)
                result_pd.append((n+1,len(rs_greed[0][0][1].split('_')),'True Optimal'))
        result_pd=pd.DataFrame(result_pd,columns=['n','Num. overlays used','method'])
    return result_pd

def get_parser():
    """Get parser object for script calculate_population_coverage.py."""

    parser = ArgumentParser()
    req_argument = parser.add_argument_group('required arguments')
    parser.add_argument("-c","--cutoff",dest='coverage_cutoff', type=float, default=6,
                        help="Target coverage lower bond when stop beam search")
    parser.add_argument("-r","--maxround",dest='max_round', type=int, default=15,
                        help="Max number of peptides to include")
    parser.add_argument("-o","--outdir", type=str, default='result',
                        help="path for results.")
    parser.add_argument("-w","--nworker", type=int, default=60,
                        help="Number of workers if use parallelization.")
    parser.add_argument("-s","--seed", type=int,
                        help="Number of workers if use parallelization.")
    return parser

if __name__ == "__main__":
    args = get_parser().parse_args()
    #result2=run_compare([args.seed],30,10,cutoff=args.cutoff,coverage=1,max_round=args.maxround,nworker=args.nworker)
    result_ex=run_compare2([args.seed],10,5,cutoff=args.cutoff,coverage=1,max_round=args.maxround,nworker=args.nworker)
    result_ex.to_pickle(args.outdir)
