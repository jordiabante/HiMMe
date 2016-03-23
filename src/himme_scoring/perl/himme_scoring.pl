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
use List::Util qw(sum reduce);
use Benchmark qw(cmpthese);
use POSIX;

# Read arguments
my $scriptname = $0;            # Get script name
my $tm_file = @ARGV[0];         # Transition matrix file name
my $ep_file = @ARGV[1];         # Emission probabilities file name
my $fasta_file = @ARGV[2];      # Get target FASTA file name
my $kmer_size= @ARGV[3];        # Get user k-mer size
my $n_entries=@ARGV[4];         # Number of entries in the FASTA file
my $outfile_sum= @ARGV[5];      # Output file

# For testing purposes
my $tm_outfile="out/tm_test.txt.gz";  # File to check the transition matrix
my $ep_outfile="out/ep_test.txt.gz";  # File to check the emission probs

# Handlers
my $FASTA;                      # Fasta file handler
my $TM;                         # Transition matrix handler
my $EP;                         # Emission probabilities handler

# Variables
my $dim=1;                      # Dimension P (based on k-mer size)
my @markov_matrix=();           # Transition Matrix [row][column]
my @emission_matrix=();         # Emission Matrix [row][column]
my $n_proc=0;                   # Number of entries processed
my $n_mem=0;                    # Number of entries stored in RAM
my $n_limit=500;                # Limit for number of entries in RAM
my $still_working=1;            # Flag
my $pi_0=1/(4**$kmer_size);     # Initial probability distribution

# Numerical stability
my @scaling_factors=();         # c_t factor in Rabiner (1989)

# Time stamps
my $st_time=0;                  # Start time
my $end_time=0;                 # End time
my $current_time=0;             # Current time
my $elapsed_time=0;             # Time elapsed

# Hashes
my %fasta_hash=();              # Hash containing sequence info of each sample
my %transition_hash=();         # Hash containing transition matrix
my %emission_hash=();           # Hash containing emission probs
my %emission_all_hash=();       # Hash containing all comb emission probs
my %score_hash_1=();            # Hash containing scores iter n-1
my %score_hash_2=();            # Hash containing scores iter n
my %results_hash=();            # Hash containing scores of all sequences
my %template_hash=();           # Hash generated based on all the combinations

################################ Main #########################################

## Read input files
read_tm();
read_ep();

## Initialize state hash and emission all hash
initialize();

## Run algorithm
run_algorithm();

# Save stuff
save_results();

############################### Subs ##########################################
## Run forward and Viterbi algorithms
sub run_algorithm
{
    # Process in chunks of n_limit
	while($still_working)
    {
        # Read n_lim entries
        read_fasta();
        # Get chromosomes from fasta
        my @entries=keys %fasta_hash;
        foreach my $entry (@entries)
        {
            # For some reason there are empty strings at the end
            next if $entry eq '';
            # Get first observation
            my $seq_1 = join '',@{$fasta_hash{$entry}}[0 .. $kmer_size-1];
            # Get all possible hidden states
            %score_hash_1=%template_hash;
            # Compute score for each hidden state
            foreach my $hidden_state (keys %score_hash_1)
            {
                $score_hash_1{$hidden_state}=$pi_0*$emission_all_hash{$hidden_state}{$seq_1};
            }
            # Compute scaling factor
            my $c_1=1/(sum values %score_hash_1);
            #print STDERR "$c_1\n";
            push @scaling_factors,$c_1;
            foreach my $hidden_state (keys %score_hash_1)
            {
                $score_hash_1{$hidden_state}*=$c_1;
                #print "$hidden_state\t$score_hash_1{$hidden_state}\n";
            }
            # Loop through the sequence
            my $seq_length=scalar @{$fasta_hash{$entry}};
            for(my $i=0;$i<=$seq_length-2*$kmer_size;$i+=$kmer_size)
            {
                # Get states in the iteration
                my $seq_1 = join '',@{$fasta_hash{$entry}}[$i .. $i+$kmer_size-1];
                my $seq_2 = join '',@{$fasta_hash{$entry}}[$i+$kmer_size .. $i+2*$kmer_size-1];
                # Get all possible hidden states i
                %score_hash_2=%template_hash;
                # Condition to leave out last incomplete k-mer
                if(length($seq_2) eq $kmer_size)
                {
                    # Consider every possible hidden state
                    foreach my $hidden_state_2 (keys %score_hash_2)
                    {
                        my $emission=$emission_all_hash{$hidden_state_2}{$seq_2};
                        my $sum=0;
                        foreach my $hidden_state_1 (keys %score_hash_1)
                        {
                            # Forward algorithm
                            my $row=$transition_hash{$hidden_state_1};
                            my $col=$transition_hash{$hidden_state_2};
                            my $score=$score_hash_1{$hidden_state_1}*
                                $markov_matrix[$row][$col]*$emission;
                            $sum+=$score;
                            #print STDERR "$hidden_state_1\t$hidden_state_2\t$score\n";
                        }
                        $score_hash_2{$hidden_state_2}=$sum;
                    }
                }
                else
                {
                    %score_hash_2=%score_hash_1;
                }
                # Scale and copy to t-1
                my $c_t=1/(sum values %score_hash_2);
                # print STDERR "$c_t\n";
                push @scaling_factors,$c_t;
                foreach my $hidden_state (keys %score_hash_2)
                {
                    $score_hash_1{$hidden_state}=$score_hash_2{$hidden_state}*$c_t;
                    #print "$hidden_state\t$score_hash_1{$hidden_state}\n";
                }
            }
            # Total score P(Y=y|HMM) = - sum_{c_t=1}^{m}(log(c_t))
            $results_hash{$entry}{length} = $seq_length;
            $results_hash{$entry}{logscore} -= log $_ for @scaling_factors;
            # Clean score hash and scoring factors
            %score_hash_1=();
            %score_hash_2=();
            @scaling_factors=();
            # Print to track progress
            print STDOUT "$entry\n";
        }
        # Clean fasta hash
        %fasta_hash=();
    }
}

## Initialize stuff
sub initialize
{
    # Get possible permutations of lenght k
    my @nucleotides=('A','C','G','T');
    my @permutations=variations_with_repetition(\@nucleotides,$kmer_size);
    # Get all possible hidden states
    foreach my $string_1 (@permutations)
    {
        my $length=$kmer_size-1;
        my $sequence_1 = join('', @{$string_1}[0..$length]);
        $template_hash{$sequence_1}=1;
        foreach my $string_2 (@permutations)
        {
            my $length=$kmer_size-1;
            my $sequence_2 = join('', @{$string_2}[0..$length]);
            my $emission=1;
            for(my $j=0;$j<$kmer_size;$j++)
            {
                my @state_1=split(//,$sequence_1);
                my @state_2=split(//,$sequence_2);
                my $row=$emission_hash{$state_1[$j]};
                my $col=$emission_hash{$state_2[$j]};
                $emission*=$emission_matrix[$row][$col];
            }
            $emission_all_hash{$sequence_1}{$sequence_2}=$emission;
        }
    }
}

## Save results
sub save_results
{
    my $sum=0;
    my $n=0;
    # Open output file
    open(OUT,">$outfile_sum") or die "Can't open file '${outfile_sum}' $!";
    # Loop through all sequences
    foreach my $key (keys %results_hash)
    {
        my $length = $results_hash{$key}{length};
        my $logscore = sprintf('%.3f', $results_hash{$key}{logscore});
        printf OUT "$key\t$length\t$logscore";
        print OUT "\n";
        $sum+=$results_hash{$key};
        $n++;
    }
    # Close handler
    close OUT;
}

## Read in fasta file
sub read_fasta
{
    my $entry;
    $n_mem=0;
    %fasta_hash=();
    while($n_mem<$n_limit)
    {
        my $line = <STDIN>;
        if($line =~ /^$/)
        {
            $still_working=0;
        }
        chomp($line);
        if( $line =~ />/)
        {   
            $entry=substr($line,1); # Get rid of leading ">" character
        }   
        else
        {   
            my @array = split //, $line;
            push @{$fasta_hash{$entry}},@array;
            $n_mem++;
        }   
	}   
}

## Print fasta
sub print_fasta
{
	# Print output
	foreach my $key (sort keys %fasta_hash)
	{
		my @entry=@{$fasta_hash{$key}};
		print ">$key\n";
		foreach my $nuc (@entry)
		{
		    print "$nuc";
		}
		print "\n";
	}
}
## Read in transition matrix
sub read_tm
{
    my @nucleotides=('A','C','G','T');
	for(my $i=1;$i<=$kmer_size;$i++)
	{
		$dim*=(scalar @nucleotides);
	}
	# Because we start with dim=0 in the loops
	$dim-=1;
    # Open file
    open($TM, "gunzip -c ${tm_file} | ")
        or die "Can't open file '${tm_file}' $!";
    my $i=0;
    # Loop through each row
    while(my $line=<$TM>)
    {
        if($i gt 0)
        {
            chomp($line);
            my @row=split(/\t/,$line);
            my $j=0;
            foreach my $col (@row)
            {
                if($j eq 0)
                {
                    $transition_hash{$col}=$i-1;
                }
                else
                {
                    $markov_matrix[$i-1][$j-1]=$col;
                }
                $j++;
            }
        }
        $i++;
    }
    close $TM;
}

## Print transition matrix
sub print_markov_matrix
{
    open(OUT, "|gzip -c > ${tm_outfile}")            # $$ is our process id
        or die "Can't open file '${tm_outfile}' $!";
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
		printf OUT "%.5f\t", $markov_matrix[$i][$j];
	    }
	    print OUT "\n";
	}
    close OUT;
}

## Read in emission probabilities
sub read_ep
{
    # Open file
    open($EP, "gunzip -c ${ep_file} | ")
        or die "Can't open file '${ep_file}' $!";
    my $i=0;
    # Loop through each row
    while(my $line=<$EP>)
    {
        if($i gt 0)
        {
            chomp($line);
            my @row=split(/\t/,$line);
            my $j=0;
            foreach my $col (@row)
            {
                if($j eq 0)
                {
                    $emission_hash{$col}=$i-1;
                }
                else
                {
                    $emission_matrix[$i-1][$j-1]=$col;
                }
                $j++;
            }
        }
        $i++;
    }
    close $EP;
}

## Print emission probabilities
sub print_emission_matrix
{
    open(OUT, "|gzip -c > ${ep_outfile}")            # $$ is our process id
        or die "Can't open file '${ep_outfile}' $!";
	# Reverse hash
	my %reverse_hash = reverse %emission_hash;
    # Print to OUT
	print OUT "\t";
	for(my $j=0;$j<=3;$j++)
	{
	    my $key = $reverse_hash{$j};
	    print OUT "$key\t";
	}
	print OUT "\n";
	for(my $i=0;$i<=3;$i++)
	{   
	    my $key = $reverse_hash{$i};
	    print OUT "$key\t";
	    for(my $j=0;$j<=3;$j++)
	    {
		printf OUT "%.5f\t", $emission_matrix[$i][$j];
	    }
	    print OUT "\n";
	}
    close OUT;
}

## log10 function
sub log10 
{
    my $n = shift;
    return log($n)/log(10);
}
##############################################################################
