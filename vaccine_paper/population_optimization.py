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

class PopulationCoverage(object):
    def __init__(self,input_epitope,allele_map,country_list,outdir):
        self.input_epitope_affinity=input_epitope
        self.allele_map=allele_map
        self.set_country(country_list)
        self.outdir=outdir
    
    def beam_search(self,beam_size=20,cutoff=0.9,max_round=1000,roll=True,curr_beam={}):
        current_max_coverage=curr_beam.items()[0][1] if len(curr_beam)>0 else 0.0
        beams=[]
        curr_length=len(curr_beam.items()[0][0].split('_')) if len(curr_beam)>0 else 0
        while current_max_coverage<cutoff and curr_length<max_round:
            cnt=curr_length
            print 'beamsearch round ',curr_length
            next_beam={}
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x.split('_') for x in curr_beam.keys()]
            pbar = tqdm(total=len(curr_candidates)*len(self.input_epitope_affinity.index))
            progress=0
            for curr_epitopes in curr_candidates:
                for next_epitope in self.input_epitope_affinity.index:
                    if roll:
                        pbar.update(1)
                    else:
                        progress+=1
                        if progress%5000==0:
                            print('scanned '+str(progress))
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
            save_pickle(join(self.outdir,'beam_'+str(cnt)+'.p'),curr_beam.items())
            save_beams(join(self.outdir,'beam_'+str(cnt)),curr_beam)
            print 'current beam: ',curr_beam.items()
            curr_length=len(curr_beam.items()[0][0].split('_'))

        print 'Coverage cutoff reached, final solution:',curr_beam.items()[0]
        print 'Per region details:'
        details=self.overall_coverage(epitopes=curr_beam.items()[0][0].split('_'))
        return curr_beam.items()[0],details,beams

    def beam_search_parallel(self,beam_size=20,cutoff=0.9,max_round=20,curr_beam={},nworker=20,diverse_cut=3):
        print 'Using %d workers' % nworker
        current_max_coverage=curr_beam.items()[0][1] if len(curr_beam)>0 else 0.0
        beams=[]
        curr_length=len(curr_beam.items()[0][0].split('_')) if len(curr_beam)>0 else 0
        pool=pp._ProcessPool(processes=nworker)
        while current_max_coverage<cutoff and curr_length<max_round:
            cnt=curr_length
            print 'beamsearch round ',curr_length
            t0 = time.time()
            self.next_beam=KthLargest(k=beam_size)
            if not len(curr_beam):
                curr_candidates=[[]]
            else:
                curr_candidates=[x.split('_') for x in curr_beam.keys()]
            print 'current candidates:',curr_candidates
            test_epitopes=[x+[y] for x in curr_candidates for y in self.input_epitope_affinity.index if diff1d(x,y,cut=diverse_cut)]
            print len(test_epitopes)
            r=pool.map_async(self.batch_scan,test_epitopes,callback=self.batch_update)
            r.get()
            curr_beam=self.next_beam.return_dict()
            current_max_coverage=curr_beam.items()[0][1]
            beams.append(curr_beam)
            save_pickle(join(self.outdir,'beam_'+str(cnt)+'.p'),curr_beam.items())
            save_beams(join(self.outdir,'beam_'+str(cnt)),curr_beam)
            print 'current beam: ',curr_beam.items()
            t1 = time.time()
            print('time passed: {}'.format(t1-t0))
            curr_length=len(curr_beam.items()[0][0].split('_'))
        pool.close()
        print 'Coverage cutoff reached, final solution:',curr_beam.items()[0]
        print 'Per region details:'
        details=self.overall_coverage(epitopes=curr_beam.items()[0][0].split('_'))
        return curr_beam.items()[0],details,beams

    def batch_scan(self,test_set):
        test_set=sorted(test_set)
        key='_'.join(test_set)
        coverage=self.overall_coverage(test_set,verbose=False)
        return coverage,key

    def batch_update(self,results):
        print 'called update'
        print results[0]
        for item in results:
            self.next_beam.add(item)

    def overall_coverage(self,epitopes,verbose=True):
        country_coverage=[]
        country_detail={}
        if isinstance(epitopes,str):
            epitopes=[epitopes]
        if len(epitopes)==1:
            binding=self.input_epitope_affinity.loc[epitopes[0]]
        else:
            single_binding=self.input_epitope_affinity.loc[epitopes]
            binding=1-(1-single_binding).product(axis=0)
        for country in self.country_list:
            if verbose:
                print country
            result=self.compute_coverage(country,binding,verbose=verbose)
            country_coverage.append(result[0])
            if verbose:
                print 'overall coverage:',country_coverage[-1]
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
            print 'locus binding probability:',locus_binding
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

def get_parser():
    """Get parser object for script calculate_population_coverage.py."""

    parser = ArgumentParser()

    req_argument = parser.add_argument_group('required arguments')
    parser.add_argument("-m","--method", type=str, default='puffin',
                        help="which prediction model to use")
    parser.add_argument("-p","--prediction", type=str, default='pred_affinity',
                        help="which type of output metric to use")
    parser.add_argument("-t","--type", type=str, default='mhc1',
                        help="which mhc type to use")
    parser.add_argument("-b","--binarize",dest='binary_cutoff', type=float,
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
    parser.add_argument("-r","--maxround",dest='max_round', type=int, default=15,
                        help="Max number of peptides to include")
    parser.add_argument("-f","--frequency",dest='freq_file', type=str, default='iedb',
                        help="which frequency file to use")
    parser.add_argument("-o","--outdir", type=str, default='result',
                        help="path for results.")
    parser.add_argument("-w","--nworker", type=int,
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
            f.writelines('%s\t%f\t%f\t%f\t%f\n' % (country,max_detail[country][0],max_detail[country][1],max_detail[country][2],max_detail[country][3]))

def save_beams(fname,beam):
    with open(fname,'w') as f:
        for seq,val in beam.items():
            f.writelines('%s\t%f\n' % (seq,val))

def save_pickle(fname,data):
    with open(fname,'w') as f:
        cPickle.dump(data,f)

def load_pickle(fname):
    with open(fname,'rb') as f:
        data=cPickle.load(f)
    return data

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
    d1=pd.DataFrame(country_detail).transpose().reset_index()
    d1.columns=['Region','Overall coverage','loci-A','loci-B','loci-C']
    d1.loc[15]=d1.mean(axis=0)
    d1.loc[15,'Region']='Average'
    ax=sns.barplot(x='Region',y='Overall coverage',data=d1,ax=axes)
    ax.set_ylim(0,1)
    for item in ax.get_xticklabels():
        item.set_rotation(90)
        item.set_fontsize(10)
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
    cv,det=pc.overall_coverage(epitopes=epitopes)
    plot_country(det,name+' peptide set',a0)
    #print cv
    model.loc[epitopes].T.reset_index(level=1, drop=True).fillna(0).plot(kind='bar',stacked=True,\
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
    cv,det=pc.overall_coverage(epitopes=epitopes)
    plot_country(det,name+' peptide set {:.2f}% coverage'.format(cv*100),a0)
    detail=model.loc[epitopes].copy()
    pos=objectives['protein'].loc[epitopes].values
    rename=['{}({})'.format(epitopes[i],pos[i]) for i in range(len(pos))]
    detail.index=rename
    detail.T.reset_index(level=1, drop=True).fillna(0).plot(kind='bar',stacked=True,\
                               colormap=ListedColormap(sns.color_palette(pal, len(epitopes))),ax=a1)    
    lg=a1.legend(bbox_to_anchor=(1, 1),ncol=len(epitopes)//12+1)
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

    if args.type=='mhc1':
        if args.freq_file=='iedb':
            frequency=pd.read_pickle('IEDB_population_frequency2392_normalized.pkl')
            #frequency=pd.read_pickle('IEDB_population_frequency102_normalized.pkl')
        else:
            print("frequency file type unknown")
            sys.exit(0)
    elif args.type=='mhc2':
        if args.freq_file=='iedb':
            #frequency=pd.read_pickle('IEDB_population_frequency_mhc2_normalized.pkl')
            frequency=pd.read_pickle('IEDB_population_frequency_mhc2_275normalized.pkl')
        else:
            print("frequency file type unknown")
            sys.exit(0)
    else:
        print("prediction output type unknown")
        sys.exit(0)

    if not exists(args.regions):
        print("Region file not exist: %s" % args.regions)
        sys.exit(0)
    regions=np.loadtxt(args.regions,dtype='S',delimiter='\n')

    if args.method not in ['puffin','deepligand','mhcflurry','netmhc','netmhc2','test','mean','mean2']:
        print("prediction method unknown")
        sys.exit(0)
    if args.prediction not in ['pred_prob','ic50','pred_affinity']:
        print("prediction output type unknown")
        sys.exit(0)

    pred_file='_'.join([args.method, args.type, args.prediction, 'pivot.pkl'])
    if not exists(pred_file):
        print("prediction file not exist %s" % pred_file)
        sys.exit(0)
    print('reading prediction file: %s' % pred_file)
    prediction=pd.read_pickle(pred_file)
    prediction_original=prediction.copy()
    if args.binary_cutoff:
        print("Binarizing predictions with threshold %f" % args.binary_cutoff)
        if args.binary_cutoff>1:
            prediction=(prediction<args.binary_cutoff).applymap(lambda x:int(x))
        else:
            prediction=(prediction>args.binary_cutoff).applymap(lambda x:int(x))
        idx = pd.IndexSlice
        prediction.loc[:,idx[:,'unknown']]=0.0
    elif args.truncate_cutoff:
        print("Truncating predictions with threshold %f" % args.truncate_cutoff)
        prediction=prediction.applymap(lambda x:0 if x<args.truncate_cutoff else x)
    
    if prediction.max().max()>1 or prediction.min().min()<0:
        print("prediction values exceeded normal range (%f,%f)" % (prediction.min().min(),prediction.max().max()))
        sys.exit(0)
    
    if len(prediction.columns)>len(frequency.columns):
        prediction=prediction[frequency.columns.values]
        prediction.columns=frequency.columns
    elif len(prediction.columns)<len(frequency.columns):
        print("mismatch between allele numbers (%d,%d)" % (len(prediction.columns),len(frequency.columns)))
        sys.exit(0)
    print("loaded %d alleles" % (len(frequency.columns)-3))

    objectives=pd.read_pickle('AllEpitopeFeatures.pkl')
    selfp=np.loadtxt('self_pept.txt',dtype='S')
    prediction.drop(selfp,errors='ignore',inplace=True)
    prediction=prediction[objectives.loc[prediction.index]['glyco_probs']<=args.glyco_cutoff]
    if args.protein:
        prediction=prediction[objectives.loc[prediction.index]['protein'].apply(lambda x:(args.protein in x))]
    prediction=prediction[objectives.loc[prediction.index]['crosses_cleavage']==0]
    prediction=prediction[objectives.loc[prediction.index]['perc_mutated']<=args.mutation_cutoff]
    prediction=prediction[prediction.sum(axis=1)>0]

    if not exists(args.outdir):
        makedirs(args.outdir)
        curr_beam={}
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
            else:
                if exists(join(args.outdir,'best_result.txt')):
                    score=float(np.loadtxt(join(args.outdir,'best_result.txt'),dtype='S')[1])
                    if score>=args.coverage_cutoff:
                        print("Optimization already reached")
                        sys.exit(0)
                    tmstamp=datetime.now().strftime("%Y%m%d%H%M%S")
                    rename(join(args.outdir,'best_result.txt'),join(args.outdir,'best_result.txt.bk'+tmstamp))
                    if exists(join(args.outdir,'plots')):
                        rename(join(args.outdir,'plots'),join(args.outdir,'plots_bk_')+tmstamp)
                beamf=join(args.outdir,'beam_{}.p'.format(max(bnum)))
                curr_beam=OrderedDict(cPickle.load(open(beamf, "rb" )))
        else:
            curr_beam={}

    if not exists(join(args.outdir,'plots')):
        makedirs(join(args.outdir,'plots'))

    pc=PopulationCoverage(prediction,frequency,regions,args.outdir)
    if not args.restart:
        prediction.to_pickle(join(args.outdir,'processed_prediction.pkl'))
        max_result=pc.overall_coverage(epitopes=list(prediction.index))
        save_detail(join(args.outdir,'upper_bound.txt'),max_result,has_overall=True)
    if args.nworker:
        best_solution,detail,beam_history=pc.beam_search_parallel(beam_size=args.beam_size,cutoff=args.coverage_cutoff,max_round=args.max_round,\
        curr_beam=curr_beam,nworker=args.nworker,diverse_cut=args.diversity)
    else:
        best_solution,detail,beam_history=pc.beam_search(beam_size=args.beam_size,cutoff=args.coverage_cutoff,max_round=args.max_round,roll=args.unroll,curr_beam=curr_beam)
    with open(join(args.outdir,'best_result.txt'),'w') as f:
        f.writelines('%s\t%f' % best_solution)
    save_detail(join(args.outdir,'best_detail.txt'),detail)

    print('plotting for peptide set:')
    print(best_solution[0].split('_'))
    plot_detail(prediction,best_solution[0].split('_'),pc,join(args.outdir,'plots'))
    plot_detail(prediction_original,best_solution[0].split('_'),pc,join(args.outdir,'plots'),'_raw',norm=(args.prediction!='ic50'))
    plot_detail3(prediction, best_solution[0].split('_'), pc,join(args.outdir,'plots'),lgd=len(best_solution[0].split('_'))<20)
    plot_detail3(prediction_original, best_solution[0].split('_'), pc, join(args.outdir,'plots'),suffix='_raw',lgd=len(best_solution[0].split('_'))<20)
    if best_solution[1]<args.coverage_cutoff:
        np.savetxt(join(args.outdir,'failed'),[best_solution[1],args.coverage_cutoff],fmt='%s')
