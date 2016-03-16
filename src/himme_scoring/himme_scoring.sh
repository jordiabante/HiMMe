#!/usr/bin/env bash
# ------------------------------------------------------------------------------
##The MIT License (MIT)
##
##Copyright (c) 2016 Jordi Abante
##
##Permission is hereby granted, free of charge, to any person obtaining a copy
##of this software and associated documentation files (the "Software"), to deal
##in the Software without restriction, including without limitation the rights
##to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
##copies of the Software, and to permit persons to whom the Software is
##furnished to do so, subject to the following conditions:
##
##The above copyright notice and this permission notice shall be included in all
##copies or substantial portions of the Software.
##
##THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
##IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
##FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
##AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
##LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
##OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
##SOFTWARE.
# ------------------------------------------------------------------------------
shopt -s extglob

abspath_script="$(readlink -f -e "$0")"
script_absdir="$(dirname "$abspath_script")"
script_name="$(basename "$0" .sh)"

# Find perl scripts
perl_script="${script_absdir}/perl/${script_name}.pl"

if [ $# -eq 0 ]
    then
        cat "$script_absdir/${script_name}_help.txt"
        exit 1
fi

TEMP=$(getopt -o hd:t:k: -l help,outdir:,threads:,kmer_size: -n "$script_name.sh" -- "$@")

if [ $? -ne 0 ] 
then
  echo "Terminating..." >&2
  exit -1
fi

eval set -- "$TEMP"

# Defaults
outdir="$PWD"
threads=2
kmer_size=1

# Options
while true
do
  case "$1" in
    -h|--help)
      cat "$script_absdir"/${script_name}_help.txt
      exit
      ;;  
    -d|--outdir)
      outdir="$2"
      shift 2
      ;;  
    -t|--threads)
      threads="$2"
      shift 2
      ;;  
    -k|--kmer_size)
      kmer_size="$2"
      shift 2
      ;;  
    --) 
      shift
      break
      ;;  
    *)  
      echo "$script_name.sh:Internal error!"
      exit -1
      ;;  
  esac
done

# Print LICENSE
cat "${script_absdir}/../../LICENSE"

# Start time
SECONDS=0
echo "[$(date)]: Starting HiMMe..."

# Inputs
tm_file="$1"
ep_file="$2"
fasta_file="$3"

# Output
fasta_basename="$(basename "$fasta_file")"
prefix="${fasta_basename%%.*}"
outprefix="${outdir}/${prefix}"
outfile="${outdir}/${prefix}_himme${kmer_size}.txt"

# Output directory
mkdir -p "$outdir"

# Count number of entries in FASTA
echo "[$(date)]: Counting number of contigs..."
n_entries="$(zcat -f "$fasta_file" | grep "^>" | wc -l)"
n_per_file="$(( $n_entries / $threads))"
echo "[$(date)]: Total of "$n_entries" contigs..."

# Split input in number of threads
echo "[$(date)]: Parallelizing..."
split_fasta.sh -d "$outdir" -n "$threads" -- "$fasta_file"
fasta_files="$(ls "$outprefix"*.fa)"

# Time elapsed
time_elapsed="$SECONDS"
echo "[$(date)]: Running HiMMe..."

# Export variables
export perl_script
export tm_file
export ep_file
export fasta_file
export kmer_size
export n_per_file
export outfile
export HIMME_PROC=0

# Run
echo "$fasta_files" | xargs -I {file} --max-proc "$threads" bash -c  \
    'zcat -f '{file}' | '$perl_script' '$tm_file' '$ep_file' '$fasta_file' '$kmer_size' '$n_per_file' '{file}.himme.tmp'' 2>&1 \
    | pv -l -s "$n_entries" -w 90 --force --timer --progress --eta --rate --interval 1 > /dev/null

# Time elapsed
time_elapsed="$SECONDS"

# Remove FASTA temporary files
echo "[$(date)]: Collapsing output HiMMe..."
rm "${outprefix}"_*.fa

# Merge output files
cat "$outprefix"*tmp | grep -v '^[[:space:]]' | sort -rg -k 2,2 > "$outfile"

# Remove temporary files
echo "[$(date)]: Cleaning..."
rm "${outprefix}"*.himme.tmp

# Time elapsed
time_elapsed="$SECONDS"
echo "[$(date)]: Total time elapsed: $(( $time_elapsed / 3600)) h $(( ($time_elapsed / 60) % 60)) m $(( $time_elapsed % 60 )) s."
