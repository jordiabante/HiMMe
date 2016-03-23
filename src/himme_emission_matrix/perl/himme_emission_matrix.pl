#!/usr/bin/env perl
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

# Libraries
use strict;
use Algorithm::Combinatorics qw(combinations variations_with_repetition);

# Read arguments
my $scriptname = $0;            # Get script name
my $vcf_file = @ARGV[0];        # Get target FASTA file name
my $n_entries=@ARGV[1];         # Number of entries in the FASTA file
my $n_bases=@ARGV[2];           # Number of bases
my $outfile=@ARGV[3];           # Output file

# Variables
my $VCF;                        # Fasta file handler
my $dim=1;                      # Dimension P (based on k-mer size)
my @emission_matrix=();         # 3D array [assembly][row][column]
my $n_proc=0;                   # Number of entries processed
my $n_mem=0;                    # Number of entries stored in RAM
my $n_limit=2000;               # Limit for number of entries in RAM
my $kmer_size=1;                


# Time stamps
my $st_time=0;                  # Start time
my $end_time=0;                 # End time
my $current_time=0;             # Current time
my $elapsed_time=0;             # Time elapsed

# Hashes
my %vcf_hash=();                # Hash containing sequence info of each sample
my %transition_hash=();         # Hash containing duplets combinations

################################ Main #########################################

# Read in fasta file
$st_time = localtime;
print STDERR "${st_time}: VCF file: ${vcf_file}\n";
# Initialize matrix
$current_time = localtime;
print STDERR "${current_time}: Initializing emission matrix ...\n";
initialize();
# Count frequencies
$current_time = localtime;
print STDERR "${current_time}: Learning emission matrix $n_bases...\n";
fill_emission_matrix();
# Print stuff
$current_time = localtime;
print STDERR "${current_time}: Saving emission matrix in ${outfile}...\n";
print_emission_matrices();

############################### Subs ##########################################
## Fill the emission matrix
sub fill_emission_matrix
{
    # Process in chunks of n_limit
	while($n_proc<$n_entries)
    {
        # Read n_lim entries
        read_vcf();
        # Get new entries
        my @entries=keys %vcf_hash;
        # Loop through the entries
        foreach my $entry (@entries)
        {
            # Get variants
            my @ref=@{$vcf_hash{$entry}{'ref'}}; 
            my @alt=@{$vcf_hash{$entry}{'alt'}}; 
            # Get index form emission_matrix
            my $row=$transition_hash{@ref[0]};
            my $col=$transition_hash{@alt[0]};
            # Fill emission_matrix
            $emission_matrix[$row][$col]+=1; 
            # Get rid of that entry 
            delete $vcf_hash{$entry};
            # Update progress
            $n_mem--;
            $n_proc++;
            my $perc=$n_proc/$n_entries*100;
            printf STDERR "\rCurrent progress: %.2f%", $perc;
        }
    }
    printf STDERR "\n";
    # Add e(i|i) and scale matrix
    for(my $i=0;$i<=$dim;$i++)
    {
        my $sum=0;
        $emission_matrix[$i][$i]=$n_bases/4;
        # Count frequencies per row
        for(my $j=0;$j<=$dim;$j++)
        {
            $sum+=$emission_matrix[$i][$j];
        }
        # Normalize each row
        if($sum!=0)
        {
            for(my $j=0;$j<=$dim;$j++)
            {
                $emission_matrix[$i][$j]/=$sum;
            }
        }
    }
}

## Initialize emission_matrices
sub initialize
{
	# Nucleotides taken into consideration
	my @nucleotides=('A','C','G','T');
	for(my $i=1;$i<=$kmer_size;$i++)
	{
		$dim*=(scalar @nucleotides);
	}
	# Because we start with dim=0 in the loops
	$dim-=1;
	# Get all possible combinations of kmer_size nucleotides
	my @permutations=variations_with_repetition(\@nucleotides,$kmer_size);
	# Codify numerically each possible permutation
	my $i=0;
	foreach my $combination (@permutations)
	{   
		my $length=$kmer_size-1;
		my $sequence = join('', @{$combination}[0..$length]);
		$transition_hash{$sequence}=$i;
		$i++;
	} 
	# Initialize emission_matrix
	for(my $i=0;$i<=$dim;$i++)
	{
	    for(my $j=0;$j<=$dim;$j++)
		{
			$emission_matrix[$i][$j]=0;
		}
	}
}

## Read in vcf file
sub read_vcf
{
    my $entry;
	while(($n_mem<$n_limit) and (my $line = <STDIN>))
	{   
		chomp($line);
		if( $line =~ /#/)
		{   
            # nothing
		}   
		else
		{   
		    my @array = split('\t',$line);
            my $id = $array[2];
            my @ref = split(//,$array[3]);
            my @alt = split(//,$array[4]);
            if((scalar @ref eq 1) and (scalar @alt eq 1))
            {
                $vcf_hash{$id}{'ref'}=[@ref];
                $vcf_hash{$id}{'alt'}=[@alt];
                $n_mem++;
            } 
            else
            {
                $n_entries--;
            }
		}   
	}   
}

## Print emission_matrices
sub print_emission_matrices
{
    open(OUT, "|gzip -c > ${outfile}")            # $$ is our process id
        or die "Can't open file '${outfile}' $!";
	# Reverse hash
	my %reverse_hash = reverse %transition_hash;
    # Print to OUT
	print OUT "\t";
	for(my $j=0;$j<=$dim;$j++)
	{
	    my $key = $reverse_hash{$j};
	    print OUT "$key\t";
	}
	print OUT "\n";
	for(my $i=0;$i<=$dim;$i++)
	{   
	    my $key = $reverse_hash{$i};
	    print OUT "$key\t";
	    for(my $j=0;$j<=$dim;$j++)
	    {
		printf OUT "%.5f\t", $emission_matrix[$i][$j];
	    }
	    print OUT "\n";
	}
    close OUT;
}
##############################################################################
