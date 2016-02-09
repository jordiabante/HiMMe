#/usr/bin/env bash

rm -rf out

set -x
../himme_ep.sh -d out -b 12 -- input.vcf.gz
