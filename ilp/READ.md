Running the ILP algorithm on a vaccine design is optimized for use through the Terminal (or a bash file).  The file containing the code is the ilp_executable.py. However, the data needs to be preprocessed with code contained within the ilp_preprocessing.ipynb notebook. The purpose of this notebook is to transform a file containing a table of peptide vs allele hits to a file of peptide vs haplotype hits. Once the new file has been created all that remains is writing the correct bash script and running it. 

Example bash script:

python ilp_executable.py -o result -fname mhc1_n3_mp30 -freq haplotype_frequency_marry.pkl -over overlap_haplotypes_mhc1.pkl -n 3 -maxpep 30 -cut 3 -maxrun 10800 -maxgap 0.1


Explanation of this bash script: 
-o is the output directory to which the results are saved.
-fname is the name of the .txt file to which the results are saved. Hence, the results are saved within the result/mhc1_n3_mp30.txt file in the above example.
-freq is the file which contains the haplotype frequencies for White, Asians, and Black. The haplotype frequencies are then averaged within the ilp_executable.py file.
-over is the overlap file that contains a table of peptide vs haplotype hits which was computed using the ilp_preprocessing.ipynb notebook.
-n is the hyperparameter n used in the maximum n-times coverage objective function.
-maxpep is the maximum number of peptides allowed in a vaccine design.
-cut If it is greater than 0, it means that the only peptides that can be contained within a vaccine design must be sufficiently different. The larger the value, the more different they must be.
-maxrun is the maximum number of seconds we allow the ILP algorithm to run for once the ILP problem has been set up by the solver. The Python-MIP solver is used. It is a solver hyperparameter.
-maxgap is the size of the step that the solver takes when solving an ILP. It is a solver hyperparameter.