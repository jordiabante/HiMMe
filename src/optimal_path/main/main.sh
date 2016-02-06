#/usr/bin/env bash

rm -rf out

set -x
../scoring_fasta.sh -k 1 -d out \
    grch37_cdna_tm2.txt.gz all_ep2991694175.txt.gz input.fa
../scoring_fasta.sh -k 2 -d out \
    grch37_cdna_tm2.txt.gz all_ep2991694175.txt.gz input.fa
../scoring_fasta.sh -k 3 -d out \
    grch37_cdna_tm3.txt.gz all_ep2991694175.txt.gz input.fa
