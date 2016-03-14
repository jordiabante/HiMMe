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
my $n_limit=100;               # Limit for number of entries in RAM

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
	while($n_proc<$n_entries)
    {
        # Read n_lim entries
        read_fasta();
        # Get chromosomes from fasta
        my @entries=keys %fasta_hash;
        foreach my $entry (@entries)
        {
            my $n_iter=floor((scalar @{$fasta_hash{$entry}})/$kmer_size)-1;
            my $pos=1;
            # Loop through the sequence
            for(my $i=0;$i<=$n_iter*$kmer_size;$i+=$kmer_size)
            {
                # Get states in the iteration
                my $seq_1=@{$fasta_hash{$entry}}[$i];
                my $seq_2=@{$fasta_hash{$entry}}[$i+$kmer_size];
                for(my $j=1;$j<$kmer_size;$j++)
                {
                    $seq_1=$seq_1.@{$fasta_hash{$entry}}[$i+$j];
                    $seq_2=$seq_2.@{$fasta_hash{$entry}}[$i+$kmer_size+$j];
                }
                # Get all possible hidden states
                %score_hash_2=%template_hash;
                # In case it's the first iteration, initialize alpha
                if($i eq 0)
                {
                    # Get all possible hidden states
                    %score_hash_1=%template_hash;
                    # Compute score for each hidden state
                    foreach my $hidden_state (keys %score_hash_1)
                    {
                        my $pi_0=1/(4**$kmer_size);
                        my $emission=$emission_all_hash{$hidden_state}{$seq_1};
                        my $score=$pi_0*$emission;
                        $score_hash_1{$hidden_state}=$score;
                    }
                    $pos++;
                }
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
                        }
                        $score_hash_2{$hidden_state_2}=$sum;
                    }
                    $pos++;
                }
                else
                {
                    %score_hash_2=%score_hash_1;
                }
                %score_hash_1=%score_hash_2;
            }
            # Total score P(X=x|HMM)
            my $total_score=0;
            foreach my $kmer (keys %score_hash_1)
            {
                $total_score+=$score_hash_1{$kmer};
            }
            # Store in results hash
            $results_hash{$entry}=$total_score;
            # CLean score hash
            %score_hash_1=();
            %score_hash_2=();
            # Update progress
            $n_proc++;
            print STDERR "$entry\n";
            #$ENV{'HIMME_PROC'}++;
            #progress_bar($ENV{'HIMME_PROC'});
            #my $perc=$ENV{'HIMME_PROC'}/$n_entries*100;
            #printf STDERR "\rCurrent progress: %.2f%", $perc;
        }
    }
}

## Initialize stuff
sub initialize
{
    # Get possible permutations of lenght k
    my @nucleotides=('A','C','T','G');
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
    foreach my $key (sort keys %results_hash)
    {
        my $string = sprintf('%.3e', $results_hash{$key});
        printf OUT "$key\t$string";
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
    my @nucleotides=('A','C','T','G');
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

# Simple progress bar
sub progress_bar
{
    my $n=shift;
    my $progress=int($n*100/$n_entries);
    my $done="";
    my $left="";
    for(my $i=1;$i<=$progress;$i++)
    {
        $done.="*";
    }
    for(my $i=$progress+1;$i<=100;$i++)
    {
        $left.="-";
    }
    print STDERR "\r$progress % [${done}${left}]";
}

##############################################################################
