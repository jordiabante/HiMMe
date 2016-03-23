#/usr/bin/env bash

rm -rf out

set -x
../run_himme.sh -k 5 -d out \
    grch37_cdna_tm5.txt.gz emission_test.txt.gz input.fa
