See Google doc here for up to date instructions: https://docs.google.com/document/d/1_RuQZ5X0tdSBYJVMp7QmrGaQY4QNC2MXyGctpXdHI-A/edit?usp=sharing

First, we need to create our protein alignments for each protein in SARS-CoV-2.  Here are the steps to do so:

- Go to the GISAID website and sign up (https://www.gisaid.org/).  Once you are approved and given a username and password, log in and go to the EpiCoV tab.  Click on Downloads and then download the FASTA and metadata files under the Genomic Epidemiology section.  Unzip these (use something like gzip or gunzip in the Terminal or Command Line) and save them to a preferred folder.  These are all the SARS-CoV-2 sequences sequenced so far.

- Familiarize yourself with Nextstrain.  Go to their SARS-CoV-2 tutorial (https://docs.nextstrain.org/en/latest/tutorials/SARS-CoV-2/steps/setup.html) and do the perform the steps below.

- Run the command: <<   git clone https://github.com/nextstrain/ncov.git    >>   in your Terminal to download the ncov folder to your local directory.

- Follow the rest of the instructions there, so run the commands << curl http://data.nextstrain.org/nextstrain.yml --compressed -o nextstrain.yml    >>
<<   conda env create -f nextstrain.yml   >>
<< conda activate nextstrain >>
<< npm install --global auspice >>
This will create a nextstrain env in your anaconda (or anaconda3).  
Open the ncov folder << cd ncov >>.
Unzip the directory example_sequences.fasta << gzip -d -c data/example_sequences.fasta.gz > data/example_sequences.fasta >> (not really necesssary for us)

- Go to the my_profiles folder within ncov.  Create a new profile (let's call it new_attempt for now).  Copy the builds.yaml and config.yaml files from the getting_started profile.  In the builds.yaml file, you can optionally change the name of the build (let's call the build my_build_attempt for now).  The default should name the build global, so that's the line you can change.  In the config.yaml file under the configfile properties change the path to the new build within your new_attempt profile you have created.  Also, under the config properties change the paths to the sequences and metadata, so that they point to the unzipped data you downloaded from GISAID.  Depending on the strength of the machine you are running the build on, you will probably want to increase the maximum allowable cores to 16 or even more.  So in the config.yaml file change cores:4 to cores:16 or perhaps even more.

- If you run the build now, there is a hugh probability that you will get an error in the augur refine job after the augur tree job has been run, especially when the dataset you are working with is large.  To avoid this, open workflow/snakemake_rules/main_workflow.smk, go to the << augur refine >> job and change the << --coalescent {params.coalescent} >> to << --coalescent opt >>.

- Create a bash file with the command << snakemake --cores 16 --profile ./my_profiles/new_attempt >> (without the << >>).  Change permissions to this file through chmod or similar commands to allow everyone to read and execute it.  Then, run the bash file on the cluster (or locally).

WARNING: This takes a significant amount of time to run.  For ~500,000 sequences it took me around 2.5 days of straight runtime.  If the process fails at some point, don't worry too much, as progress is periodically saved.  Just rerun the bash file after figuring out the error.

- Now, you need to translate your genes to proteins and get protein files with aligned proteins for each protein (ex. E, M, N, S, ORF1a, ORF8, etc.)

- Create a new bash file (or run the command locally in the Terminal) and include the command: augur translate --tree results/my_build_attempt/tree.nwk --ancestral-sequences results/my_build_attempt/nt_muts.json --reference-sequence sars_cov2.gff --output-node-data my_outputs/aa_muts.json --alignment-output my_outputs/aligned_protein_%GENE.fasta

where sars_cov2.gff is the latest gff file for SARS-CoV-2 you can download from the web (ex. https://www.ncbi.nlm.nih.gov/sars-cov-2/)

- Run another augur translate process but instead of sars_cov2.gff use the default nextstrain default in defaults/reference_seq.gb within the ncov folder.  Compare the 2 results and merge because there are some proteins that are only present in one of the two augur translate runs.



 
(Alex's) Notes on the Python notebooks (ipynb files):

- ParseProteins.ipynb: (No need to run it) Not super important, it simply takes in a fasta file containing all the SARS proteins (ORF1a, M, ... etc.) and outputs a csv file (SARS_all_protein_epitopes.csv) that is a pandas dataframe of protein fragment, which proteins (ORF1a, M, ... etc.) they came from, as well as the starting position of the fragment within each protein 0-indexed (data/processed/SARS_all_protein_epitopes.csv).

- ListofCleavageSites.ipynb: (No need to run it) Outputs a small csv file (data/processed/ProteinCleavageSites.csv) containing the cleavage sites for ORF1a, ORF1b (ORF1ab is ORF1a concatenated with ORF1b), and the Spike Protein S.  The csv file contains the overall protein unit (ORF1a, ORF1b, S) where the smaller protein comes from, along with its name (ex. Guanine-N7 methyltransferase) and the start and end positions of the amino acids within the larger protein overall unit 0-indexed.

- GlycosylationPreds.ipynb: (No need to run it) Outputs a csv file (data/SARS_glycosolation_sites_processed.csv) where all the glycosylation sites are reported.  There are 90 of them according to Trenton's data and each entry corresponds to one glycosylation site and includes: the protein (ex. 1a, 3a, 9b, E, ... etc.), the position within the protein that is glycosylated, the sequence that is glycosylated (4 amino acids), and the probability of glycosylation. Furthermore, it does some Glycosylation site processing.

- MapEpitopestoProteins:  (No need to run it) Nothing important, related to vaccine epitopes.

IMPORTANT:

- CurrentMutationRate.ipynb: 1st step - run this first with your new protein alignments

- EpitopeMutationRate.ipynb: 2nd step - run this next

- CombiningAllAdditionalVaccineConsiderations: 3rd step - run this afterwards


