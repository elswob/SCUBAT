#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

#7 args
#1. -t target fasta file
#2. -q query fasta file
#3. -p psl BLAT output
#4. -o overlap length (max 10 due to possibility of small blat alignments allowing incorrect associations)
#5. -l length of alignment
#6. -c number of processors
#7. -m mismatch score

my %EST_hash=();
my $target;
my $query;
my $start;
my $stop;

my $target_file = "";
my $query_file = "";
my $psl = "";
my $overlap = 10;
my $cores = 1;
my $cov_length = 90;
my $joiner = 'cap3';
my $mismatch = 90;

GetOptions (
	"target=s" => \$target_file,
	"query=s" => \$query_file,
	"psl=s" => \$psl,
	"overlap=s" =>\$overlap,
	"cores=s" =>\$cores,
	"length=s" =>\$cov_length,
	"mismatch=s" =>\$mismatch,
);

if ($target_file eq "" || $query_file eq "" || $psl eq "" || $overlap >20){
	die("SCUBAT.pl \n-t contigs file \n-q transcripts file \n-p blat psl output \n-o overlap length (default 10, max 20)\n-c number of processors (default 1)\n-l min \% coverage length of transcript (default 90)\n-m mismatch score (default 0.92)\n");
} 


my $dir = "o.$overlap.lc.$cov_length.j.$joiner.m.$mismatch";
unless (-d "$dir"){
	`mkdir $dir`;
}


#adjust the mismatch value
$mismatch = $mismatch / 100;

open S,">$dir/run_log.txt";
my $cov_length_div = $cov_length / 100;

####################################################################################################################################################

print "Getting target lengths...\n";
print S "Getting target lengths...\n";
open C, $target_file or die;
my $header=""; my $len = "";
my %target_length;
my %target_seq;
while(<C>){
	chomp;
	#forward slashes and pipes in seq names were giving problems later on
	$_ =~ s/\//_/g;
	$_ =~ s/\|/_/g;
	if ($_ =~ m/^>.*/){
		if ($header ne ""){
			$target_length{$header}=length($len);
			$target_seq{$header}=$len;
		}	
		if ($_ =~ m/^>(.*?)\s+.*/){
			$header = $1;
		}elsif ($_ =~ m/^>(.*)/){
			$header = $1;
		}
		$len = "";
	}else{
		$len = $len.$_;
	}
}
#catch the last one
$target_length{$header}=length($len);
$target_seq{$header}=$len;

open TS,">$dir/target_seqs.fa";
for (sort keys %target_seq){
	print TS "$_\n$target_seq{$_}\n";
}
#open TD,">$dir/target_dump.txt";
#print TD Dumper( \%target_seq );

####################################################################################################################################################

print "Getting query lengths...\n";
print S "Getting query lengths...\n";
open Q, $query_file or die;
$header=""; 
$len = "";
my %query_length;
my %query_seq;
while(<Q>){
	chomp;
	#forward slashes and pipes in seq names were giving problems later on
	$_ =~ s/\//_/g;	
	$_ =~ s/\|/_/g;
	if ($_ =~ m/^>(.*?)/){
		if ($header ne ""){
			$query_length{$header}=length($len);
			$query_seq{$header}=$len;
		}
		if ($_ =~ m/^>(.*?)\s+.*/){
			$header = $1;
		}elsif($_ =~ m/^>(.*)/){
			$header = $1;
		}
		$len = "";
	}else{
		$len = $len.$_;
	}
}
#catch the last one
$query_length{$header}=length($len);
$query_seq{$header}=$len;

open QS,">$dir/query_seqs.fa";
for (sort keys %query_seq){
	print QS "$_\n$query_seq{$_}\n";
}
#open QD,">$dir/query_seq.txt";
#print QD Dumper( \%query_seq );

####################################################################################################################################################

print "Reading blat output...\n";
print S "Reading blat output...\n";
open MIS,">$dir/mismatch.txt";
my $num_blat = `grep -c '' $psl`;
chomp $num_blat;
my $count_blat=0;
open B, $psl or die;
while (<B>){
	$count_blat++;
	my $counter = $count_blat / $num_blat * 100;
	my $rounded = sprintf("%.2f", $counter);
	unless ($counter % 1){
		printf("\r%9d percent processed",$rounded);	
	}
	#forward slashes and pipes in seq names were giving problems later on
	$_ =~ s/\//_/g;
	$_ =~ s/\|/_/g;
	my @array = split("\t");
	#check if blat results contain more contigs than given in file
	next unless exists ($target_seq{$array[13]});
	#check matches - mismatches / matches is greater than something
	if (($array[0] - $array[1]) / $array[0] < $mismatch){
		print MIS $_;
		next;
	}
	my $or = $array[8];
	$target = $array[13];
	unless (exists $target_seq{$target}){
		die("A target seq in the BLAT output is not in the target file!\nPerhaps you have the wrong orientation?\n");
	}
	$query = $array[9];
	$start = $array[11];
	$stop = $array[12];
	
	#sort the start and stop positions on the subject
	my @pos = ($start,$stop);
	my @sort_pos = sort {$a <=> $b} @pos;
	push(@sort_pos, $or); 
	#check if contig ID is already in hash array
	#if more than one hit per query then take the two extreme positions
	#create a hash of hashes with the sorted positions as values in an array 
	if (exists $EST_hash{$query}{$target}){
		if ($start < $EST_hash{$query}{$target}[0]){
			$EST_hash{$query}{$target}[0] = $start;
		}if ($stop > $EST_hash{$query}{$target}[1]){
			$EST_hash{$query}{$target}[1] = $stop;
		}
	}else{
		$EST_hash{$query}{$target} = [@sort_pos];
	}
}
#open EH,">$dir/EST_hash.txt";
#print EH Dumper( \%EST_hash );

####################################################################################################################################################

#create groups of contigs per EST
#only include if there are no ambiguous contigs 
print "\nMaking initial groups...\n";
print S "Making initial groups...\n";
my %groups;
my %one_wrong_groups;
my %over_one_wrong_groups;
my $contigID;

unless (-d "$dir/one_contig/"){
	`mkdir $dir/one_contig`;
}
my $count_single=0;
open OC,">$dir/one_contig/single_map.fa";
my $count_est = 0;
foreach my $EST_value (sort keys %EST_hash){
	$count_est++;
	my $counter = $count_est / (keys %EST_hash) * 100;
	my $rounded = sprintf("%.2f", $counter);
	unless ($counter % 1){
		printf("\r%9d percent processed",$rounded);	
	}
	breakpoint: while($EST_value){
		if (keys (%{$EST_hash{$EST_value}}) > 1){
			my $stop = 0;
			my $check = 0;
			#print "--- $EST_value --- \n";
			for $contigID (sort {$EST_hash{$EST_value}{$a}[0]<=>$EST_hash{$EST_value}{$b}[0]} keys %{$EST_hash{$EST_value}}){
				#print "contigID = $contigID - start = $EST_hash{$EST_value}{$contigID}[0] - stop = $EST_hash{$EST_value}{$contigID}[1]\n";
				if ($EST_hash{$EST_value}{$contigID}[0] >= $stop-$overlap){
					$groups{$EST_value}{$contigID} = [@{$EST_hash{$EST_value}{$contigID}}];
					$stop = $EST_hash{$EST_value}{$contigID}[1];
				}else{
					#allow one non-consecutive contig as it may be allelic contigs
					delete $groups{$EST_value};
					$stop = 0;
					for $contigID (sort {$EST_hash{$EST_value}{$a}[0]<=>$EST_hash{$EST_value}{$b}[0]} keys %{$EST_hash{$EST_value}}){
						if ($EST_hash{$EST_value}{$contigID}[0] >= $stop-$overlap){
							$one_wrong_groups{$EST_value}{$contigID} = [@{$EST_hash{$EST_value}{$contigID}}];
							$stop = $EST_hash{$EST_value}{$contigID}[1];		
						}else{
							$one_wrong_groups{$EST_value}{$contigID} = [@{$EST_hash{$EST_value}{$contigID}}];
							$check++;
						}
					}
					if ($check > 1){
						delete $one_wrong_groups{$EST_value};
						for $contigID (sort {$EST_hash{$EST_value}{$a}[0]<=>$EST_hash{$EST_value}{$b}[0]} keys %{$EST_hash{$EST_value}}){
							$over_one_wrong_groups{$EST_value}{$contigID} = [@{$EST_hash{$EST_value}{$contigID}}];
						}
						last breakpoint;
					}else{
						last breakpoint;
					}
					
				}	
			}
			last breakpoint;
		}
		#catch the single contigs that capture the whole transcript
		elsif(keys (%{$EST_hash{$EST_value}}) == 1){
			for my $contigID (keys %{$EST_hash{$EST_value}}){
				next unless exists $query_length{$EST_value};
				#print "$EST_hash{$EST_value}{$contigID}[1] - $EST_hash{$EST_value}{$contigID}[0] / $query_length{$EST_value} > $cov_length_div\n";
				if ((($EST_hash{$EST_value}{$contigID}[1] - $EST_hash{$EST_value}{$contigID}[0]) / $query_length{$EST_value}) > $cov_length_div){
					print OC ">$contigID.$EST_value\n$target_seq{$contigID}\n";
					$count_single++;
				}
			}
			last breakpoint;
		}else{
			last breakpoint;
		}
	} 
}

#empty the EST_hash
undef %EST_hash;

print "\n$count_single transcripts mapped completely to one contig covering over $cov_length percent of the transcript.\n"; 
print S "\n$count_single transcripts mapped to one contig over $cov_length percent of the transcript.\n"; 
print "Created ".(keys %groups)." groups.\n";
print S "Created ".(keys %groups)." groups.\n";
#check if 'enough' of the transcript is covered
my $count_short=0;
for my $k1 (keys %groups){
	my $min=100000000;
	my $max=0;
	#get start and finish coordinates and check they are grater than x % of transcript
	for $contigID (sort {$groups{$k1}{$a}[0]<=>$groups{$k1}{$b}[0]} keys %{$groups{$k1}}){
		if ($groups{$k1}{$contigID}[0] < $min){
			$min = $groups{$k1}{$contigID}[0];
		}if ($groups{$k1}{$contigID}[1] > $max){
			$max = $groups{$k1}{$contigID}[1];
		}
	}
	if (exists $query_length{$k1}){
		my $c = ($max - $min / $query_length{$k1});
		if ((($max - $min) / $query_length{$k1}) > $cov_length_div){
			next;
		}else{
			$count_short++;
			#print "$k1 does not cover enough of the transcript!\n";
			delete($groups{$k1});
		}
	}
}
print "$count_short groups were removed for covering less than $cov_length percent of the transcript.\n";
print S "$count_short groups were removed for covering less than $cov_length percent of the transcript.\n";

#create group id headers to include contig IDs
my %groups_header;
for my $k1 (keys %groups){
	my $cat_contig="";
	for my $k2 (keys %{$groups{$k1}}){
		$cat_contig.=$k2.$groups{$k1}{$k2}[2];
	}
	$groups_header{$k1}=$cat_contig;
}

#print "Counting contigs in good groups...\n";
my %group_contigs;
for my $k1 (keys %groups){
	#print "$k1\n";
	for my $k2 (%{$groups{$k1}}){
		if (exists $groups{$k1}{$k2}){
			#print "$k2\n";
			$group_contigs{$k2}++;
		}
	}
}
print "There are ".(keys %groups)." groups with no problems and ".(keys %group_contigs)." contigs involved\n";
print S "There are ".(keys %groups)." groups with no problems  and ".(keys %group_contigs)." contigs involved\n";

open ONE,">$dir/one_wrong_groups_$overlap.txt";
for my $k1 (sort keys %one_wrong_groups){
	next unless exists $query_length{$k1};
	print ONE "--- $k1 [$query_length{$k1}] ---\n";
	for my $k2 (sort {$one_wrong_groups{$k1}{$a}[0]<=>$one_wrong_groups{$k1}{$b}[0]} keys %{$one_wrong_groups{$k1}}){
		print ONE "$k2 [$target_length{$k2}]->$one_wrong_groups{$k1}{$k2}[0]\t$one_wrong_groups{$k1}{$k2}[1]\t$one_wrong_groups{$k1}{$k2}[2]\n";
	}
}

#print "Counting contigs in one wrong groups...\n";
my %one_wrong_group_contigs;
for my $k1 (keys %one_wrong_groups){
	#print "$k1\n";
	for my $k2 (%{$one_wrong_groups{$k1}}){
		if (exists $one_wrong_groups{$k1}{$k2}){
			#print "$k2\n";
			$one_wrong_group_contigs{$k2}++;
		}
	}
}

#print "Counting contigs in over one wrong groups...\n";
my %over_one_wrong_group_contigs;
for my $k1 (keys %over_one_wrong_groups){
	#print "$k1\n";
	for my $k2 (%{$over_one_wrong_groups{$k1}}){
		if (exists $over_one_wrong_groups{$k1}{$k2}){
			#print "$k2\n";
			$over_one_wrong_group_contigs{$k2}++;
		}
	}
}

print "There are ".(keys %one_wrong_groups)." groups of contigs with one error with ".(keys %one_wrong_group_contigs)." contigs involved\n";
print S "There are ".(keys %one_wrong_groups)." groups of contigs with one error with ".(keys %one_wrong_group_contigs)." contigs involved\n";
print "There are ".(keys %over_one_wrong_groups)." groups of contigs with more than one error with ".(keys %over_one_wrong_group_contigs)." contigs involved\n";
print S "There are ".(keys %over_one_wrong_groups)." groups of contigs with more than one error with ".(keys %over_one_wrong_group_contigs)." contigs involved\n";

####################################################################################################################################################

$count_short=0;
print "Checking groups with one error...\n";
print S "Checking groups with one error...\n";
for my $k1 (keys %one_wrong_groups){
	next unless exists $query_length{$k1};
	my $min=100000000;
	my $max=0;
	#don't want less than two
	if (keys %{$one_wrong_groups{$k1}} > 2){
		#get start and finish coordinates and check they are grater than x % of transcript
		for $contigID (sort {$one_wrong_groups{$k1}{$a}[0]<=>$one_wrong_groups{$k1}{$b}[0]} keys %{$one_wrong_groups{$k1}}){
			if ($one_wrong_groups{$k1}{$contigID}[0] < $min){
				$min = $one_wrong_groups{$k1}{$contigID}[0];
			}if ($one_wrong_groups{$k1}{$contigID}[1] > $max){
				$max = $one_wrong_groups{$k1}{$contigID}[1];
			}
		}
		if ((($max - $min) / $query_length{$k1}) > $cov_length_div){
			my $stop=0;
			my @prev;
			#print "--- $k1 ---\n";
			#go back through and pick the best contig	
			for $contigID (sort {$one_wrong_groups{$k1}{$a}[0]<=>$one_wrong_groups{$k1}{$b}[0]} keys %{$one_wrong_groups{$k1}}){
				if ($one_wrong_groups{$k1}{$contigID}[0] >= $stop-$overlap){
					#print "good $k1 - $contigID - $one_wrong_groups{$k1}{$contigID}[0] $one_wrong_groups{$k1}{$contigID}[1]\n";
					$groups{$k1}{$contigID} = [@{$one_wrong_groups{$k1}{$contigID}}];
					$stop = $one_wrong_groups{$k1}{$contigID}[1];
					@prev=($contigID,$one_wrong_groups{$k1}{$contigID}[0],$one_wrong_groups{$k1}{$contigID}[1]);
				}else{
					#if the conflicting contig covers more of the transcript use that one
					#print "bad $k1 - $contigID - $one_wrong_groups{$k1}{$contigID}[0] $one_wrong_groups{$k1}{$contigID}[1]\n";
					if (($one_wrong_groups{$k1}{$contigID}[1] - $one_wrong_groups{$k1}{$contigID}[0]) > ($prev[2] - $prev[1])){
						#print "replacing $prev[0] - $prev[1] - $prev[2] with... \n$k1 - $contigID - $one_wrong_groups{$k1}{$contigID}[0] $one_wrong_groups{$k1}{$contigID}[1]\n";
						delete $groups{$k1}{$prev[0]};
						$groups{$k1}{$contigID} = [@{$one_wrong_groups{$k1}{$contigID}}];
					}elsif (($one_wrong_groups{$k1}{$contigID}[0] == $prev[1]) && ($one_wrong_groups{$k1}{$contigID}[1] == $prev[2])){
						#print "same as previous, keep longest contig...\n";
						#print "prev $prev[0] = $target_length{$prev[0]}\n new $contigID = $target_length{$contigID}\n";
						if ($target_length{$contigID} > $target_length{$prev[0]}){
							delete $groups{$k1}{$prev[0]};
							$groups{$k1}{$contigID} = [@{$one_wrong_groups{$k1}{$contigID}}];
							#print "use new\n";
						}else{
							#print "use previous\n";
						}
					}else{
						#print "keeping $prev[0] - $prev[1] - $prev[2]\n";
					}
				}
			}
						
		}else{
			$count_short++;
		}
	}
}

#print "Counting contigs in good groups...\n";
%group_contigs=();
for my $k1 (keys %groups){
	#print "$k1\n";
	for my $k2 (%{$groups{$k1}}){
		if (exists $groups{$k1}{$k2}){
			#print "$k2\n";
			$group_contigs{$k2}++;
		}
	}
}

print "$count_short groups with one wrong contig were removed for covering less than $cov_length percent of the transcript.\n";
print S "$count_short groups with one wrong contig were removed for covering less than $cov_length percent of the transcript.\n";
print "There are now ".(keys %groups)." good groups with ".(keys %group_contigs)." contigs involved\n";
print S "There are now ".(keys %groups)." good groups with ".(keys %group_contigs)." contigs involved\n";
#open GR,">$dir/groups_$overlap.txt";
#print GR Dumper( \%groups );

####################################################################################################################################################

#open OVER,">$dir/over_one_wrong_groups_$overlap.txt";
#print OVER Dumper( \%over_one_wrong_groups );

open OUT,">$dir/scaffold_run_$overlap.txt";

unless (-d "$dir/merged_contigs"){
	`mkdir $dir/merged_contigs`;
}
print "Merging ".(keys %groups)." groups based on shared contigs...\n";
print S "Merging ".(keys %groups)." groups based on shared contigs...\n";
my %merge_table;
my $count_groups=0;
my $id_count=1;
my $m_count=1;
my %join_groups; my %join_groups_final; my $cat_keys; my $group_flag; my $join_groups_count=0;
for my $key1 (sort keys %groups){
	if (keys %{$groups{$key1}} > 1){
		my $tmp_merge = "";
		$count_groups++;
		my $counter = $count_groups / (keys %groups) * 100;
		my $rounded = sprintf("%.2f", $counter);
		unless ($counter % 1){
			printf("\r%9d percent groups processed",$rounded);	
		}
		#print "------- $key1 [$query_length{$key1}]--------\n";
		#check query seqs in blat exist in query file
		next unless exists $query_seq{$key1};
		print OUT "------- $key1 [$query_length{$key1}]--------\n";
		#open the new merged files
		open MERGE,">$dir/merged_contigs/$key1.merge.fa";
		$tmp_merge .= ">scaffold_m"."$m_count\n";
		$merge_table{$key1}="scaffold_m".$m_count;
		$m_count++;
		for my $key2 (sort {$groups{$key1}{$a}[0]<=>$groups{$key1}{$b}[0]} keys %{$groups{$key1}}){
			$group_flag = "";
			#add to join_groups
			for my $k1 (keys %join_groups){ 
				$cat_keys = "";
				if (exists $join_groups{$k1}{$key2}){
					for my $key22 (sort {$groups{$key1}{$a}[0]<=>$groups{$key1}{$b}[0]} keys %{$groups{$key1}}){
						$join_groups{$k1}{$key22}++;
						$cat_keys.= $key22."-";  
					}
					$cat_keys =~ s/(.*?)-$/$1/;
					$join_groups_final{$k1}{$key1} = $cat_keys;
					$group_flag = "yes";
				}
			}
			print OUT "$key2 [$target_length{$key2}]-> ";
			foreach (@{$groups{$key1}{$key2}}){
				print OUT "$_\t";
			}
			print OUT "\n";
			#print "key2 = $key2\n";
			#print new merged contig files and rev comp where necessary (use regex as blat can sometimes use ++ and +-)
			if (@{$groups{$key1}{$key2}}[2] !~ /-/){
				$tmp_merge .= "$target_seq{$key2}NNNNNNNNNN";
			}else{
				#revcomp those with a -
				my $revcomp = reverse($target_seq{$key2});
				$revcomp =~ tr/ACGTacgt/TGCAtgca/;
				$tmp_merge .= "$revcomp"."NNNNNNNNNN";
			}	
		}
		if ($group_flag ne "yes"){
			$join_groups_count++;
			for my $key22 (sort {$groups{$key1}{$a}[0]<=>$groups{$key1}{$b}[0]} keys %{$groups{$key1}}){
				$join_groups{$join_groups_count}{$key22}++;
				$cat_keys.= $key22."-";  
			}
			$cat_keys =~ s/(.*?)-$/$1/;
			$join_groups_final{$join_groups_count}{$key1} = $cat_keys;
		}

		$id_count++;
		#remove the last set of NNNNNNNNNNNs
		$tmp_merge =~ s/(.*?)NNNNNNNNNN$/$1/;
		print MERGE "$tmp_merge\n";
		close MERGE;
	}
}

#no longer need %groups, %target_length, %query_length or %query_seq; 
undef %groups;
undef %target_length;
undef %query_length;
undef %query_seq;

#open JG,">$dir/join_groups.txt";
#print JG Dumper( \%join_groups);

#open PRE_JGF,">$dir/pre_join_groups_final_$overlap.txt";
#print PRE_JGF Dumper( \%join_groups_final );

#join_groups_final hash can be redundant due to a contig group sharing contigs with separate groups
#need to merge these groups together

my (%bin,$current_bin,%finalgroup,%match,$match_bin);

open RED,">$dir/red_check.txt";
for my $id ( sort keys %join_groups_final ) {
    print RED "--- $id ---\n";
    $current_bin++;
    #need to make hash of overlapping group ids as each group could contain ids from separate groups
    #use this list to then make one group containing all the linked transcripts and contigs
    %match = ();
    for my $trans (keys %{ $join_groups_final{$id} }) {
    	print RED "transcript = $trans\n";
        if ( exists $bin{$trans} ) {
        	$match{$bin{$trans}}=""; 
        	print RED "<- transcript $trans exists ->\n";
        	#$current_bin = $bin{$trans};
        	$match_bin = $bin{$trans};
        	print RED "cur bin = $current_bin\n";
        	print RED "match bin = $match_bin\n";
        	#last;
        }
    }if ((keys %match) > 0){
    	    print RED "*** Making new group ***\n";
    	    for my $join_id (keys %match){  	
    	    	    for my $trans (keys %{ $finalgroup{$join_id} }) {
    	    	    	    print RED "$match_bin - $trans = join_id: $join_id - $trans\n";
    	    	    	    $finalgroup{$match_bin}{$trans} = $finalgroup{$join_id}{$trans};
    	    	    	    $bin{$trans} = $match_bin;
    	    	    }
    	    	    unless ($join_id eq $match_bin){
    	    	    	    delete($finalgroup{$join_id});
    	    	    }
    	    }
    	    for my $trans (keys %{ $join_groups_final{$id} }) {
    	    	    print RED "$match_bin - $trans = id: $id - $trans\n";
    	    	    $finalgroup{$match_bin}{$trans} = $join_groups_final{$id}{$trans};
    	    	    $bin{$trans} = $match_bin;
    	    }
    }else{
    	    for my $trans (keys %{ $join_groups_final{$id} }) {
    	    	    $bin{$trans} = $current_bin;
    	    	    print RED "new bin = $bin{$trans} - $trans\n";
    	    	    $finalgroup{$current_bin}{$trans} = $join_groups_final{$id}{$trans};
    	    }
    }
}

#just to print it all out
for $current_bin (keys %finalgroup) {
    print RED "--- ".$current_bin ." ---\n";
    for my $trans (keys %{$finalgroup{$current_bin}}){
    	    print RED "$trans - $finalgroup{$current_bin}{$trans}\n";
    }
}


%join_groups_final = %finalgroup;


print "\nBased on shared contigs the groups form ".(keys %join_groups_final)." super groups\n";
print S "Based on shared contigs the groups form ".(keys %join_groups_final)." super groups\n";

#open JGF,">$dir/join_groups_final_$overlap.txt";
#print JGF Dumper( \%join_groups_final );


####################################################################################################################################################

#=pod

#join up the groups
unless (-d "$dir/join"){
	`mkdir $dir/join`;
}
print "Creating directories and files for each merged group...\n";
print S "Creating directories and files for each merged group...\n";
my $caps;
my $count_caps=0;
my $pwd = `pwd`;
open JOIN,">$dir/join_commands.txt";
chomp $pwd;
unless (-d "$dir/singles"){
	`mkdir $dir/singles`;
}
for my $k1 (sort {$a<=>$b} keys %join_groups_final){
	$count_caps++;
	my $counter = $count_caps / (keys %join_groups_final) * 100;
	my $rounded = sprintf("%.2f", $counter);
	unless ($counter % 1){
		printf("\r%9d percent joined groups processed",$rounded);	
	}
	#print "$k1\n";
	if (keys %{$join_groups_final{$k1}} > 1){
		unless (-d "$dir/join/$k1"){
			#`mkdir $dir/join/$k1`;
			system("mkdir $dir/join/$k1");
		}
		for my $k2 (keys %{$join_groups_final{$k1}}){
			#`ln -s $pwd/$dir/merged_contigs/$join_groups_final{$k1}{$k2}.merge.fa $dir/join/$k1`;
			unless (-e "$dir/join/$k1/$k2.merge.fa") {
				system("ln -s $pwd/$dir/merged_contigs/$k2.merge.fa $dir/join/$k1");
			}
		}
		#print "cat $dir/join/$k1/*.fa > $dir/join/$k1/$k1.cat\n";
		#`cat $dir/join/$k1/*.fa > $dir/join/$k1/$k1.cat`;
		system("cat $dir/join/$k1/*.fa > $dir/join/$k1/$k1.cat");
		#print CAP "cat $dir/join/$k1/*.fa > $dir/join/$k1/$k1.cat\n";
		
		print JOIN "$joiner $dir/join/$k1/$k1.cat > $dir/join/$k1/$k1.log 2> $dir/join/$k1/$k1.err\n";
		#print CAP "phrap $dir/join/$k1/$k1.cat > $dir/join/$k1/$k1.log\n";
		
		
		#print CAP "rm $dir/join/$k1/$k1.cat\n";
	}else{
		for my $k2 (keys %{$join_groups_final{$k1}}){
			#`ln -s $pwd/$dir/merged_contigs/$join_groups_final{$k1}{$k2}.merge.fa $dir/singles`;
			unless (-e "$dir/singles/$k2.merge.fa") {
				system("ln -s $pwd/$dir/merged_contigs/$k2.merge.fa $dir/singles");
			}
		}
	}
}
print "\nRunning $joiner on merged groups using $cores processors...\n";
`parallel -P $cores < $dir/join_commands.txt`;



####################################################################################################################################################

print "Finished assembling\n";

print "Checking through $joiner joins...\n";
my %contig_fail;
my %contig_join;
my $contig_join_header;
my %contig_merge;
my $contig_split;
my %largest_merge_contig;
$count_caps=0;
if ($joiner eq 'cap3'){
	open CAPCHECK,">$dir/cap3_check.txt";
	open JF,">$dir/joins_table.txt";
	open JG,">$dir/joined_groups.fa";
	open JFX,">$dir/joins_table_failed.txt";
	for my $k1 (sort {$a<=>$b} keys %join_groups_final){
		$count_caps++;
		my $counter = $count_caps / (keys %join_groups_final) * 100;
		my $rounded = sprintf("%.2f", $counter);
		unless ($counter % 1){
			printf("\r%9d percent joined groups checked",$rounded);	
		}
		if ((keys %{$join_groups_final{$k1}}) > 1){
			my $find_largest_merge="";
			my $find_largest_merge_count=0;
			if (-s "$dir/join/$k1/$k1.cat.cap.singlets"){
				#find the merge with the most contigs
				print CAPCHECK "$k1 failed\n";
				for my $k2 (keys %{$join_groups_final{$k1}}){
					my @s = split(/-/,$join_groups_final{$k1}{$k2});
					my $merge_count=0;
					foreach (@s){
						$merge_count++;
					}
					if ($merge_count > $find_largest_merge_count){
						$find_largest_merge=$k2;
						$find_largest_merge_count = $merge_count;
					}
				}
				my %keep_contigs=();
				print CAPCHECK "keeping $find_largest_merge $merge_table{$find_largest_merge}\n";
				#create list of contigs involved in merge that is being kept
				my @s = split(/-/,$join_groups_final{$k1}{$find_largest_merge});
				foreach (@s){
					$keep_contigs{$_}="";
					$largest_merge_contig{$_}="";
				}
				for my $k2 (keys %{$join_groups_final{$k1}}){
					if ($k2 eq $find_largest_merge){
						print JF "$merge_table{$k2}\tsingle\t".$k2."\t".$join_groups_final{$k1}{$k2}."\n";
					}else{
						if (exists $merge_table{$k2}){
							#remove merged files except longest one 
							print JFX "$merge_table{$k2}\tsingle\t".$k2."\t".$join_groups_final{$k1}{$k2}."\n";
							print CAPCHECK "removing $k2 - $merge_table{$k2}\n";
							`rm $dir/merged_contigs/$k2.merge.fa`;
							delete($merge_table{$k2});
							my @s = split(/-/,$join_groups_final{$k1}{$k2});
							foreach (@s){
								unless (exists $keep_contigs{$_}){
									$contig_fail{$_}++;
								}
							}
						}
					}
				}
			}else{
				$contig_join_header="";
				for my $k2 (keys %{$join_groups_final{$k1}}){
					if (exists $merge_table{$k2}){
						#remove merged files
						print CAPCHECK "$k1 assembled -> removing $k2 - $merge_table{$k2}\n";
						`rm $dir/merged_contigs/$k2.merge.fa`;
						delete($merge_table{$k2});
						my @s = split(/-/,$join_groups_final{$k1}{$k2});
						foreach (@s){
							$contig_join{$_}++;
						}
						print JF "scaffold_j"."$k1\tgroup\t".$k2."\t$join_groups_final{$k1}{$k2}\n";
						$contig_join_header .= $k2."[$join_groups_final{$k1}{$k2}]";
					}
					
				}
				open J, "$dir/join/$k1/$k1.cat.cap.contigs";
				$contig_split=0;
				while (<J>){
					if ($_ =~ /^>.*/){
						$contig_split++;
						print JG ">scaffold_j"."$k1"."_"."$contig_split\n";
					}else{
						print JG "$_";
					}
				}
			}
		}else{
			for my $k2 (keys %{$join_groups_final{$k1}}){
				print JF "$merge_table{$k2}\tsingle\t".$k2."\t".$join_groups_final{$k1}{$k2}."\n";
				my @s = split(/-/,$join_groups_final{$k1}{$k2});
				foreach (@s){
					#print "removing $_\n";
					$contig_merge{$_}++;
				}
			}
		}
			
	}	
}


print "\n".(keys %contig_fail)." contigs failed and will be added to the final contigs file\n";
print S "\n".(keys %contig_fail)." contigs failed and will be added to the final contigs file\n";
print "".(keys %contig_join)." contigs were successfully joined as part of a group\n";
print S "".(keys %contig_join)." contigs were successfully joined as part of a group\n";
print "".(keys %contig_merge)." contigs were successfully joined as single merge events\n";
print S "".(keys %contig_merge)." contigs were successfully joined as single merge events\n";


print "Printing new contig set...\n";
open C,">$dir/tmp.contigs.fa";
for (sort keys %target_seq){
	#check if the contig was part of a single merge, a succesful join or the lonfest merge of a failed join
	unless ((exists $contig_join{$_}) || (exists $contig_merge{$_}) || (exists $largest_merge_contig{$_})){
		print C ">$_\n$target_seq{$_}\n";
	}
}
#make joined_groups.fa single line fasta
my $jg_header="";
my $jg_seq="";
open F, "$dir/joined_groups.fa";
open FOUT,">$dir/joined_groups_sl.fa";
while(<F>){
	chomp;
	if ($_ =~ /^>.*/){
		if ($jg_header ne ""){
			print FOUT $jg_header."\n$jg_seq\n";
		}
		$jg_header = $_;
		$jg_seq="";
	}else{
		s/\n//;
		$jg_seq .= $_;
	}
}
print FOUT $jg_header."\n$jg_seq\n";

`mv $dir/joined_groups_sl.fa $dir/joined_groups.fa`;
`cat $dir/tmp.contigs.fa $dir/merged_contigs/* $dir/joined_groups.fa > $dir/contigs.fa`;
`rm $dir/tmp.contigs.fa`;


