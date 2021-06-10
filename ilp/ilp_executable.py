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


from mip import * ### Mixed Integer Programming Library












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
                        help="File to read the peptide vs alleles data")
    parser.add_argument("-n", "--ntarget", type=int, default=5,
                        help="The ntarget for max n-times coverage")
    parser.add_argument("-maxpep", "--max_number_of_pepts", type=int, default=30,
                        help="The maximum number of peptides allowed in a vaccine")
    parser.add_argument("-cut", "--cut_size", type=int, default=3,
                        help="The minimum edit distance for us to consider peptides being sufficiently distinct")
    parser.add_argument("-maxrun", "--max_solver_runtime", type=int, default=3600,
                        help="The number of seconds allotted to the solver to solve the ILP problem")
    parser.add_argument("-maxgap", "--max_solver_gap", type=float, default=0.1,
                        help="The maximum gap for the solver's steps")



    
    return parser




if __name__ == "__main__":


    


    ############################################################################################################
    ############################################################################################################
    ############################################################################################################
    ### ILP Code


    args = get_parser().parse_args()

    frequency = pd.read_pickle(args.frequency)
    print('frequency = ', frequency)

    ### The haplotype frequency file considers White, Asians, and Black:  We just average them all out here.  Depending on the target
    ### population we may want to change the simple averaging to weighted averaging
    average_frequency = pd.DataFrame(index = ['Average'], columns = frequency.columns)
    average_frequency.loc['Average',:] = frequency.sum(axis=0)/3

    print('average_frequency = ', average_frequency)

    overlap = pd.read_pickle(args.overlap)
    print('overlap = ', overlap)


    outdir = args.outdir
    print('outdir = ', outdir)

    if not exists(outdir):
        makedirs(outdir)



    ### Set these values from inputs
    max_number_of_peptides = args.max_number_of_pepts
    print('max_number_of_peptides = ', max_number_of_peptides)
    n_target = args.ntarget
    print('n_target = ', n_target)
    cut = args.cut_size
    print('cut = ', cut)
    max_solver_runtime = args.max_solver_runtime
    print('max_solver_runtime = ', max_solver_runtime)
    max_gap = args.max_solver_gap
    print('max_gap = ', max_gap)



    ### https://docs.python-mip.com/en/latest/quickstart.html for a starting guide for the mip package
    m = Model()

    m = Model(sense=MAXIMIZE, solver_name=CBC)

    num_of_peptides = len(overlap.index)
    num_of_haplotypes = len(overlap.columns)

    # max_number_of_peptides = 15
    # n_target = 5


    a = [m.add_var(var_type=BINARY) for i in range(num_of_peptides)] ### The a_i's

    # print(m)

    t = [m.add_var(var_type=BINARY) for j in range(num_of_haplotypes)] ### The t_j's

    m += xsum(a[i] for i in range(num_of_peptides)) <= max_number_of_peptides, 'max_number_of_peptides'


    counter = 0
    for j in range(num_of_haplotypes):
        counter += 1
        if counter % 50 == 0:
            print('ILP Round = ', counter)
        m += -1000*(1-t[j]) <= xsum(a[i] * overlap.iloc[i,j] for i in range(num_of_peptides)) - n_target + 0.1, 'left hand side of big M for j='+str(j)
        m += xsum(a[i] * overlap.iloc[i,j] for i in range(num_of_peptides)) - n_target + 0.1 <= 1000 * t[j], 'right hand side of big M for j='+str(j)

    m.objective = maximize(xsum(t[j] * average_frequency.iloc[0,j] for j in range(num_of_haplotypes)))

    print('Model = ', m)

    print('Number of variables = ', m.num_cols)

    print('Number of Constraints = ', m.num_rows)



    # def diff1d(curr,seq,cut=3):  ### WRONG!!!
    #     for x in curr:
    #         if (x in seq) or (seq in x):
    #             if abs(len(x)-len(seq))<=cut:
    #                 return False
    #         else:
    #             for diff_l in range(1,cut):
    #                 for diff_r in range(1,cut-diff_l+1):
    #                     if x[diff_l:]==seq[:-diff_r] or seq[diff_l:]==x[:-diff_r]:
    #                         return False
    #     return True

    def diff1d(x,y,cut=3):
        if (x in y) or (y in x):
            if abs(len(x)-len(y))<=cut:
                return False
        else:
            for diff_l in range(1,cut):
                for diff_r in range(1,cut-diff_l+1):
                    if x[diff_l:]==y[:-diff_r] or y[diff_l:]==x[:-diff_r]:
                        return False
        return True



    counter = 0
    counter_redundant = 0
    for i in range(len(overlap.index)):
        for j in range(i+1, len(overlap.index)):
            counter += 1
            pept1 = overlap.index[i]
            pept2 = overlap.index[j]
            if not diff1d(pept1, pept2, cut):
                counter_redundant += 1
    #             print('Not Sufficiently different!!!')
    #             print(pept1)
    #             print(pept2)
                m += a[i] + a[j] <= 1, 'Non-redundancy condition ' + str(i) + ' ' + str(j)

    print('Number of redundant pairs = ', counter_redundant)

    print('Model after nonredundancy conditions have been added = ', m)

    print('Number of variables = ', m.num_cols)

    print('Number of Constraints = ', m.num_rows)


    output_peptide_vars = []
    output_peptide_ind = []

    counter_pept = 0

    m.max_gap = max_gap ### 0.05 or 0.1
    status = m.optimize(max_seconds=max_solver_runtime) ### I usually let it run between 30 mins (1800 seconds) to 5 hours (18000 seconds)
    # status = m.optimize(max_seconds=18000)
    if status == OptimizationStatus.OPTIMAL:
        print('optimal solution cost {} found'.format(m.objective_value))
    elif status == OptimizationStatus.FEASIBLE:
        print('sol.cost {} found, best possible: {}'.format(m.objective_value, m.objective_bound))
    elif status == OptimizationStatus.NO_SOLUTION_FOUND:
        print('no feasible solution found, lower bound is: {}'.format(m.objective_bound))
    if status == OptimizationStatus.OPTIMAL or status == OptimizationStatus.FEASIBLE:
        print('solution:')
        for v in m.vars:
            if abs(v.x) > 1e-6: # only printing non-zeros
                print('{} : {}'.format(v.name, v.x))
                if counter_pept < max_number_of_peptides:
                    temp_name = v.name
                    output_peptide_vars.append(temp_name)
                    output_peptide_ind.append(int(temp_name.split('(')[1][:-1]))
                    counter_pept += 1

    else:
        print(status)


    print('output_peptide_vars = ', output_peptide_vars)
    print('output_peptide_ind = ', output_peptide_ind)

    output_peptides = [overlap.iloc[x,:].name for x in output_peptide_ind]

    
    print('output_peptides = ', output_peptides)




    # output_peptides.to_csv(outdir + args.file_name)

    f = open(outdir + '/' + args.file_name + '.txt', 'w')
    f.write(str(output_peptides))
    f.close()














