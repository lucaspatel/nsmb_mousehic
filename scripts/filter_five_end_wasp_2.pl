#!/usr/bin/perl

#Modified version to remove secondary flag mark
#modified again by DG to preserve read order for filter_remapped_reads.py (part of WASP pipeline)

use strict;
use warnings;

my $prev_id = "";
my @five;
my @three;
my @unmap;
my @mid;
my @all;
my $counter = 0;
my $seq;
my $qual;

while (<STDIN>){
	chomp;
	if (/^@/){
		print $_."\n";
		next;
	}
	my ($id, $flag, $chr_from, $loc_from, $mapq, $cigar, $d1, $d2, $d3, $read, $read_qual, @rest) = split /\t/;
	my $bin = reverse(dec2bin($flag));
	my @binary = split(//,$bin);
	#if ($prev_id ne $id && $prev_id ne ""){				#DG EDIT HERE======================================================================
	if ($binary[11] == 0 && $prev_id ne ""){				#DG EDIT HERE======================================================================
		if ($counter == 1 && @five == 1){
			print $five[0]."\n";
		}
		elsif ($counter == 2 && @five == 1){
			my ($id_1, $flag_1, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $read_1, $read_qual_1, @rest_1) = split /\t/, $five[0];
			my $bin_1 = reverse(dec2bin($flag_1));
			my @binary_1 = split(//,$bin_1);
			$binary_1[11] = 0;
			my $bin_1_new = reverse(join("",@binary_1));
			my $flag_1_new = bin2dec($bin_1_new);
			if ($binary_1[4]==1){
				$seq = reverse ($seq);
				$seq =~tr/ATCG/TAGC/;
				$qual = reverse ($qual);
			}
			$cigar_1 =~tr/H/S/;
			print(join("\t",$id_1, $flag_1_new, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $seq, $qual, @rest_1)."\n");
		}

		#DG EDIT 1 STARTS HERE==========================================================================================================================
		else{
			my ($id_1, $flag_1, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $read_1, $read_qual_1, @rest_1) = split /\t/, $all[0];
			$chr_from_1="*";	$loc_from_1=0;
			print(join("\t",$id_1, $flag_1, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $read_1, $read_qual_1, @rest_1)."\n");
		}
		#DG EDIT 1 ENDS HERE============================================================================================================================

		$counter = 0;
		undef @unmap;
		undef @five;
		undef @three;
		undef @mid;
		undef @all;
		$seq="";
	}

	$counter++;
	$prev_id = $id;
	push @all,$_;
	if ($cigar !~ /H/){
		$seq = $read;
		$qual = $read_qual;
		if ($binary[4]==1){
			$seq = reverse ($seq);
			$seq =~tr/ATCG/TAGC/;
			$qual = reverse ($qual);
		}
	}
	if ($binary[2]==1 || $cigar =~ /[ID]/){
		push @unmap,$_;
	}
	elsif ($binary[4]==0 && $cigar =~ m/^[0-9]*M/ || $binary[4]==1 && $cigar =~ m/.*M$/){
		push @five, $_;
	}
	elsif ($binary[4]==1 && $cigar =~ m/^[0-9]*M/ || $binary[4]==0 && $cigar =~ m/.*M$/){
		push @three, $_;
	}
	elsif ($cigar =~ m/^[0-9]*[HS].*M.*[HS]$/){
		push @mid, $_;
	}
}

if ($counter == 0){
}
elsif ($counter == 1 && @five == 1){
	print $five[0]."\n";
}
elsif ($counter == 2 && @five == 1){
	my ($id_1, $flag_1, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $read_1, $read_qual_1, @rest_1) = split /\t/, $five[0];
	my $bin_1 = reverse(dec2bin($flag_1));
	my @binary_1 = split(//,$bin_1);
	$binary_1[11] = 0;
	my $bin_1_new = reverse(join("",@binary_1));
	my $flag_1_new =  bin2dec($bin_1_new);
	if ($binary_1[4]==1){
		$qual = reverse ($qual);
		$seq = reverse ($seq);
		$seq =~tr/ATCG/TAGC/;
	}
	$cigar_1 =~tr/H/S/;
	print(join("\t",$id_1, $flag_1_new, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $seq,$qual, @rest_1)."\n");
}
#DG EDIT 2 STARTS HERE==========================================================================================================================
else{
	my ($id_1, $flag_1, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $read_1, $read_qual_1, @rest_1) = split /\t/, $all[0];
	$chr_from_1="*";	$loc_from_1=0;
	print(join("\t",$id_1, $flag_1, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $read_1, $read_qual_1, @rest_1)."\n");
}
#DG EDIT 2 ENDS HERE============================================================================================================================

sub dec2bin {
	my $str = unpack("B32", pack("N", shift));
	return $str;
}

sub bin2dec {
	return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}
