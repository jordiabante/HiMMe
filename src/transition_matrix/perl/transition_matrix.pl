#!/usr/bin/env perl
# ------------------------------------------------------------------------------
##The MIT License (MIT)
##
##Copyright (c) 2015 Jordi Abante
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
use DateTime;
use Algorithm::Combinatorics qw(combinations variations_with_repetition);
use Math::Matrix;                       # Matrix operations

# Read arguments
my $scriptname = $0;            # Get script name
my $fasta_file = @ARGV[0];      # Get target FASTA file name
my $outfile = @ARGV[1];      	# Get output prefix
my $kmer_size= @ARGV[2];        # Get user k-mer size

# Variables
my $dim=1;                      # Dimension P (based on k-mer size)
my @markov_matrix=();           # 3D array [assembly][row][column]

# Time stamps
my $st_time=0;                  # Start time
my $end_time=0;                 # End time
my $current_time=0;             # Current time
my $elapsed_time=0;             # Time elapsed

# Hashes
my %fasta_hash=();              # Hash containing sequence info of each sample
my %transition_hash=();         # Hash containing duplets combinations

######################################### Main ######################################

# Read in fasta file
$st_time = DateTime->now();
print STDERR "${st_time}:Reading fasta ...\n";
read_fasta();
#print_fasta();

# Initialize Markov matrices
$current_time = DateTime->now();
print STDERR "${current_time}:Initializing Markov matrix ...\n";
initialize();
#fill_markov_matrix();

# Print stuff
#save_markov_matrices();

########################################### Subs ######################################
## Fill the markov matrix
sub fill_markov_matrix
{
	# Get chromosomes from fasta
	my @chormosomes=keys %fasta_hash;
	foreach my $chromosome (@chromosomes)
	{
		for(my $i=0;$i<=(scalar @{$fasta_hash{$chromosome}})-2*$kmer_size;$i++)
		{
		    # Get nucleotides of the iteration
		    my $seq_1=@{$fasta_hash{$chromosome}}[$i];
		    my $seq_2=@{$fasta_hash{$chromosome}}[$i+$kmer_size];
		    for(my $j=1;$j<$kmer_size;$j++)
		    {
			$seq_1=$seq_1.@{$fasta_hash{$chromosome}}[$i+$j];
			$seq_2=$seq_2.@{$fasta_hash{$chromosome}}[$i+$kmer_size+$j];
		    }
		    # Get codon index form markov_matrix
		    my $row=$transition_hash{$seq_1};
		    my $col=$transition_hash{$seq_2};
		    # Fill markov_matrix
		    $markov_matrix[$row][$col]+=1; 
		    #print STDERR "$first_codon\t$second_codon\n";
		}
		# Scale matrix
		my $sum;
		for(my $i=0;$i<=$dim;$i++)
		{
			$sum=0;
			# Count mutations per row
			for(my $j=0;$j<=$dim;$j++)
			{
				$sum+=$markov_matrix[$i][$j];
			}
			# Normalize each row
			if($sum!=0)
			{
				for(my $j=0;$j<=$dim;$j++)
				{
					$markov_matrix[$i][$j]/=$sum;
				}
			}
		}
	}
}

## Initialize markov_matrices
sub initialize
{
	# Nucleotides taken into consideration
	my @nucleotides=('A','C','T','G');
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
		# Initialize markov_matrix
		for(my $i=0;$i<=$dim;$i++)
		{
		    for(my $j=0;$j<=$dim;$j++)
		    {
			$markov_matrix[$i][$j]=0;
		    }
		}
	}
}

## Read in fasta file
sub read_fasta
{
	my $chromosome;
	open(my $FH, $fasta_file) or die "Could not open file '$fasta_file' $!";
	while( my $line = <$FH>)
	{   
		chomp($line);
		if( $line =~ />/)
		{   
		    $chromosome=substr($line,1); # Get rid of leading ">" character
		    $fasta_hash{$chromosome}="";
		}   
		else
		{   
		    my @array = split //, $line;
		    $fasta_hash{$chromosome} = [@array];
		}   
	}   
}

## Print fasta
sub print_fasta
{
	# Print output
	foreach my $key (sort keys %fasta_hash)
	{
		my @chromosome=@{$fasta_hash{$key}};
		print ">$key\n";
		foreach my $nuc (@chromosome)
		{
		    print "$nuc";
		}
		print "\n";
	}
}

## Save markov_matrices
sub save_markov_matrices
{
	open(my $fh, '>', "$outfile");
	# Reverse hash
	my %reverse_hash = reverse %transition_hash;
	print $fh "\t";
	for(my $j=0;$j<=$dim;$j++)
	{
	    my $key = $reverse_hash{$j};
	    print $fh "$key\t";
	}
	print $fh "\n";
	for(my $i=0;$i<=$dim;$i++)
	{   
	    my $key = $reverse_hash{$i};
	    print $fh "$key\t";
	    for(my $j=0;$j<=$dim;$j++)
	    {
		printf $fh "%.3f\t", $markov_matrix[$i][$j];
	    }
	    print $fh "\n";
	}
	close $fh;
}

###########################################################################################
