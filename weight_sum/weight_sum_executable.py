import numpy as np
import pandas as pd
import seaborn as sns

import itertools
import math
import scipy.special
import matplotlib.pyplot as plt
from os.path import dirname, basename,join,exists
from os import makedirs,system,listdir,rename
from copy import copy, deepcopy
from argparse import ArgumentParser


def optivax_robust(over, hap, thresh, set_of_peptides):
    '''
    Evaluates the objective value for optivax robust
    EvalVax-Robust over haplotypes!!!
    '''
    
    num_of_peptides = len(over.index)
    num_of_haplotypes = len(over.columns)
    
    
    def my_filter(my_input, thresh):
        '''
        A simple function that acts as the indicator variable in the pseudocode
        '''
        for index in range(len(my_input)):
            if my_input[index] >= thresh:
                my_input[index] = 1
            else:
                my_input[index] = 0
        return my_input



    total_overlays = np.zeros(num_of_haplotypes, dtype=int)

    for pept in set_of_peptides:

        total_overlays = total_overlays + np.array(over.loc[pept,:])

    filtered_overlays = my_filter(total_overlays, thresh)
    
#     print('filtered_overlays = ', filtered_overlays)
#     print('hap = ', hap)

    return np.sum(filtered_overlays * np.array(hap))



def diff1d(x,y,cut=3):
    '''
    A helper function that when given two peptides x and y, returns True if the peptides are sufficiently different and False
    if they are too similar. The larger the value of cut, the more different the peptides need to be in order for the function
    to return a value of True

    Inputs:
        - x: a string object that represents the one peptide using the single letter code for each amino acid
        - y: a string object that represents the other peptide using the single letter code for each amino acid
        - cut: an int object whose magnitude influences the output of the diff1d function
    Output:
        - True (sufficiently different peptides x and y) or False (not sufficiently different peptides x and y = x and y
            are too similar)
    '''
    if (x in y) or (y in x):
        if abs(len(x)-len(y))<=cut:
            return False
    else:
        for diff_l in range(1,cut):
            for diff_r in range(1,cut-diff_l+1):
                if x[diff_l:]==y[:-diff_r] or y[diff_l:]==x[:-diff_r]:
                    return False
    return True


def weight_sum_non_redundant(over, hap, max_number_of_elems, cut=3):
    '''
    Weight Sum algorithm
    Inputs:
        - over: a pandas.DataFrame object, where the rows are the peptides and the columns are the HLA haplotypes that 
            consist of 3 loci and can be MHC Type I (HLA-A, HLA-B, HLA-C) or MHC Type II (HLA-DR, HLA-DQ, HLA-DP)
        - hap: a pandas.DataFrame object, where there is only one row and the columns are the HLA haplotypes; it usually represents
            the average over (White, Asian, Black) frequency of the haplotypes; it can be changed to represent any desired
            haplotype frequency
        - max_number_of_elems: an int object, which represents the maximum number of peptides allowed in a vaccine design
        - cut: an int object (default = 3), which represents the cut value under which we consider two peptides to be too
            similar to each other to include both in a single vaccine design; the way similarity is determined is by querying 
            the function diff1d: if it returns False, it means that the two peptides are not different enough, so they are in 
            other words too similar; set equal to 0 to ignore the diff1d criterion
    Output:
        - best_subset: a list object that contains max_number_of_elems many peptides that are all sufficiently different given
            cut value. 

    '''
    sort_list = []

    for pept in over.index:
        temp_score = np.sum(np.array(over.loc[pept,:]) * np.array(hap))
        sort_list.append((pept, temp_score))
        

    sort_list.sort(key = lambda x: x[1], reverse = True)
    
    
    best_subset = []


    ### if we want to omit the step of ensuring that no two peptides in a vaccine design are too similar
    if cut == 0:
        best_subset = []
        for j in range(min(max_number_of_elems, len(over.index))):
            best_subset.append(sort_list[j][0])
        return best_subset

    
    stop = min(max_number_of_elems, len(over.index))
    counter = 0
    flag = True
    i = 0
    
    while flag:
        pept = sort_list[i][0]
        add = True
        for prev_pept in best_subset:
            if not diff1d(prev_pept, pept, cut=cut):
                add = False
                break
        if add:
            best_subset.append(sort_list[i][0])
            counter += 1
        if counter == stop:
            flag = False
            
        i +=1

    return best_subset







def get_parser():
    """Get parser object for script calculate_population_coverage.py."""

    parser = ArgumentParser()

    req_argument = parser.add_argument_group('required arguments')

    parser.add_argument("-o", "--outdir", type=str, default='result',
                        help="Path for results")
    parser.add_argument("-fname", "--file_name", type=str, default="try1",
                        help="The name the output file should have within the output directory")
    parser.add_argument("-freq", "--frequency", type=str,
                        help="File to read the haplotype frequencies from")
    parser.add_argument("-over", "--overlap", type=str,
                        help="File to read the peptide vs alleles or peptide vs haplotype data")
    parser.add_argument("-o_a", "--overlap_allele", type=int, default=0,
                        help="1 if the --overlap file passed in is peptide vs alleles and 0 if it is peptide vs haplotypes and has already been binarized")
    # parser.add_argument("-n", "--ntarget", type=int, default=5,
    #                     help="The ntarget for max n-times coverage")
    parser.add_argument("-maxpep", "--max_number_of_pepts", type=int, default=30,
                        help="The maximum number of peptides allowed in a vaccine")
    parser.add_argument("-c", "--cut", type=int, default=3,
                        help="The cut value for ommitting peptides that are too similar; a value of 0 should be provided if similar peptides are not to be excluded from a vaccine design.")


    
    return parser




if __name__ == "__main__":


    args = get_parser().parse_args()

    






    # frequency = pd.read_pickle('../optivax_saber/haplotype_frequency_marry.pkl')
    frequency = pd.read_pickle(args.frequency)
    print('frequency = ')
    print(frequency)
    
    try:
        overlap = pd.read_pickle(args.overlap)
    except:
        if args.overlap_allele == 0:
            overlap = pd.read_csv(args.overlap, index_col=0, header=[0,1,2])
        elif args.overlap_allele == 1:
            overlap = pd.read_csv(args.overlap, index_col=0, header=[0])
    print('overlap')
    print(overlap)
    # raise Exception


    if args.overlap_allele == 1:
        overlap_binarized = overlap.applymap(lambda x: 1 if x>0.638 else 0)
        zero_hit_peptides = overlap_binarized.sum(axis=1)[overlap_binarized.sum(axis=1) == 0].index ### Find the peptides with zero hits
        overlap_binarized.drop(zero_hit_peptides, axis='index', inplace=True) ### Drop the zero hit peptides from the DataFrame
        # overlap_binarized_final = overlap_binarized.droplevel('loci', axis=1) ### You might need to uncomment this and comment out the line below
        overlap_binarized_final = overlap_binarized ### Comment this out if you uncomment the above line (will depend on the format of the input data)

        ### Sum up the hits across all alleles of each haplotype
        overlap_haplotypes = pd.DataFrame(index = overlap_binarized_final.index, columns = frequency.columns)

        unique_columns = overlap_binarized_final.columns.unique()
        count = 0
        # print(unique_columns)
        for pept in overlap_haplotypes.index:

            count += 1
            if count%50 == 0:
                print('Round = ', count)

            for col in overlap_haplotypes.columns:
                temp_sum = 0
                for hap_type in col:
                    if hap_type in unique_columns:
                        temp_sum += int(overlap_binarized_final.loc[pept, hap_type])
        
                overlap_haplotypes.loc[pept, col] = temp_sum

            # if count == 100:
            #     print(overlap_haplotypes)
            #     raise Exception

        print('overlap_haplotypes')
        print(overlap_haplotypes)

    elif args.overlap_allele == 0:
        overlap_haplotypes = overlap




    ### The haplotype frequency file considers White, Asians, and Black:  We just average them all out here.  Depending on the target
    ### population we may want to change the simple averaging to weighted averaging
    average_frequency = pd.DataFrame(index = ['Average'], columns = frequency.columns)
    average_frequency.loc['Average',:] = frequency.sum(axis=0)/3

    # overlap = pd.read_csv('../optivax_saber/weight_sum_preprocessed.csv', index_col=0, header=[0,1,2])

    # outdir = 'weight_sum_outputs/try1'
    outdir = args.outdir

    # ### The value of n_target
    # # thresh = 3
    # thresh = args.ntarget

    ### The value of cut
    cut = args.cut

    ### The maximum number of peptides allowed in a vaccine
    # max_number_of_elems = 70
    max_number_of_pepts = args.max_number_of_pepts


    if not exists(outdir):
        makedirs(outdir)


    # best_subset, score = weight_sum(overlap_haplotypes, average_frequency, thresh, max_number_of_pepts)

    best_subset = weight_sum_non_redundant(overlap_haplotypes, average_frequency, max_number_of_pepts, cut)

    print('best_subset = ')
    print(best_subset)


    f = open(outdir + '/' + args.file_name + '.txt', 'w')
    f.write(str(best_subset))
    f.close()





















