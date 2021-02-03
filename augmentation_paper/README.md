# OptiVax subunit augmentation and vaccine design

This repository contains code and part of the data for the following paper:
<br>
[Predicted Cellular Immunity Population Coverage Gaps for SARS-CoV-2 Subunit Vaccines and their Augmentation by Compact Peptide Sets](https://doi.org/10.1016/j.cels.2020.11.010)
<br />
Cell Systems, 2021
<br />
Authors: Ge Liu, Brandon Carter, David K. Gifford


## Set up
Copy the `*_mhc*_pept.txt` files to the same root directory of the python files. 
Download the haplotype frequency data and model predictions for `ensemb-adapt` model from https://www.dropbox.com/sh/v1jcin4mh7jua14/AAB7W0Y7IXtXRL8Ehlrtvft6a?dl=0, 
and put them into the same root directory. Copy over the `AllEpitopeFeatures.pkl` and `self_pept.txt` from the home directory of this github repo.

## Non-redundant compression of subunit peptides
```
usage: non-redundant_compression.py [-h] [-m METHOD] [-p PREDICTION] [-t TYPE]
                                    [-b BINARY_CUTOFF] [-tr TRUNCATE_CUTOFF]
                                    [-gl GLYCO_CUTOFF] [-mt MUTATION_CUTOFF]
                                    [-s BEAM_SIZE] [-c COVERAGE_CUTOFF]
                                    [-ic INITIAL_CUT] [-lo LOWER_BOUND]
                                    [-r MAX_ROUND] [-f FREQ_FILE] [-o OUTDIR]
                                    [-w NWORKER] [-bs BATCHSIZE] [-re REGIONS]
                                    [-cr CORRECTION] [-pr PROTEIN]
                                    [-d DIVERSITY] [--unroll] [--downsamp]
                                    [--restart] [--skippre]
```
The following example compress the RBD redundant windows of MHC1 peptides into non-redundant peptdies, with the MIRA corrected ensemble model (ensemb-adapt).
```
python non-redundant_compression.py -o downsample_RBD_mhc1_beam1_1.0_normed -lo 1000 -b 0.638 \
    -t mhc1_haplotype -s 1 -c 0.99 -m ensemb-adapt -r 2000 -p pred_affinity -w 224 -bs 20 -d 3 -mt 1.0 \
    -pr RBD_mhc1_pept.txt --skippre 
```
## Augmentation with independent peptide set for MHC1 and MHC2
```
usage: augmentation_independ.py [-h] [-m METHOD] [-p PREDICTION] [-t TYPE]
                                [-b BINARY_CUTOFF] [-tr TRUNCATE_CUTOFF]
                                [-gl GLYCO_CUTOFF] [-mt MUTATION_CUTOFF]
                                [-s BEAM_SIZE] [-c COVERAGE_CUTOFF]
                                [-ic INITIAL_CUT] [-high RATIO]
                                [-low RATIO_LOW] [-lo LOWER_BOUND]
                                [-r MAX_ROUND] [-f FREQ_FILE] [-o OUTDIR]
                                [-w NWORKER] [-bs BATCHSIZE] [-re REGIONS]
                                [-pr PROTEIN] [-ba BASEMENT] [-bd BASEFILE]
                                [-pf PREDFILE] [-cr CORRECTION] [-d DIVERSITY]
                                [--unroll] [--downsamp] [--restart]
                                [--skippre]
```
The following example computes augmentation set for RBD MHC1, with the MIRA corrected ensemble model (ensemb-adapt). The result will be saved into `augment_RBD_with_all_mhc1` folder.
```
python augmentation_independ.py -o augment_RBD_with_all_mhc1 -lo 8 -b 0.638 \
    -t mhc1_haplotype -s 5 -c 0.99 -m ensemb-adapt -p pred_affinity -w 224 -bs 40 -d 3 -mt 0.001 -ba RBD_mhc1_pept.txt \
    -bd downsample_RBD_mhc1_beam1_1.0_normed_seq.txt --skippre 
```
With additional argument `-pf`, the following example computes augmentation set for RBD MHC1 using only MIRA positive candidates (or any list of candidates user specified).
```
python augmentation_independ.py -o augment_RBD_with_MIRAonly_mhc1 -lo 8 -b 0.638 \
    -t mhc1_haplotype -s 10 -c 0.99 -m ensemb-adapt -p pred_affinity -w 96 -bs 40 -d 1 -mt 1.0 -ba RBD_mhc1_pept.txt \
    -bd downsample_RBD_mhc1_beam1_1.0_normed_seq.txt -pf Adaptive_candidate_mhc1_normed38.pkl --skippre
```
