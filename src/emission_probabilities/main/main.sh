#/usr/bin/env bash

rm -rf out

set -x
../emission_probabilities.sh -d out -- input.vcf.gz
