#!/usr/bin/perl

#########################################
#Combine paired reads after wasp pipeline
#########################################

#########################################
#Modified at May 02 2016
#Change die to warn when headers are not consistent
#Implement read name comparison as samtool sort by name
#########################################

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

#CHECK INPUT_____________________________________________________________________
if(scalar(@ARGV)!=2){die "Incorrect usage;";}
my ($bam_1, $bam_2)=@ARGV;
open(FILE1,"samtools view -h $bam_1 |") or die "input file 1 error!";	#Arguement 1: Read 1 sam file
open(FILE2,"samtools view -h $bam_2 |") or die "input file 2 error!"; #Arguement 2: Read 2 sam file

my $one_only=0;
my $two_only=0;
my $both=0;
my $count=0;

my $one=<FILE1>;
my $two=<FILE2>;

while ($one =~ /^@/ && $two =~ /^@/){
	if ($one ne $two){$!=1;print STDERR "inconsistent bam header.\n";}
	else{
		print $one;
	}
	$one=<FILE1>;
	$two=<FILE2>;
	next;
}

if ($one =~ /^@/ || $two =~ /^@/){
	$! = 1;
	print STDERR "inconsistent bam header.\n";
}

my ($id_1, $flag_1, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $read_1, $read_qual_1, @rest_1) = split(/\t/,$one);
my ($id_2, $flag_2, $chr_from_2, $loc_from_2, $mapq_2, $cigar_2, $d1_2, $d2_2, $d3_2, $read_2, $read_qual_2, @rest_2) = split(/\t/,$two);

while(!(eof(FILE1) || eof(FILE2)))	
{
	$count++;
	if($count%10000000 == 0)	{
		print STDERR "### $count\n";	}

	chomp $one;
	chomp $two;


	if ($id_1 eq $id_2)	{
		combine();
		$one=<FILE1>;
		($id_1, $flag_1, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $read_1, $read_qual_1, @rest_1) = split(/\t/,$one);
		$two=<FILE2>;
		($id_2, $flag_2, $chr_from_2, $loc_from_2, $mapq_2, $cigar_2, $d1_2, $d2_2, $d3_2, $read_2, $read_qual_2, @rest_2) = split(/\t/,$two);
	}
	
	elsif (lessthan($id_1, $id_2))	{
		$one_only++;	
		$one=<FILE1>;
		($id_1, $flag_1, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $read_1, $read_qual_1, @rest_1) = split(/\t/,$one);	}
	
	else	{
		$two_only++;	
		$two=<FILE2>;
		($id_2, $flag_2, $chr_from_2, $loc_from_2, $mapq_2, $cigar_2, $d1_2, $d2_2, $d3_2, $read_2, $read_qual_2, @rest_2) = split(/\t/,$two);	}	
}

if(eof(FILE1))	{
	while(!eof(FILE2) && ($id_1 gt $id_2))	{
		if ($id_1 eq $id_2)	{
			combine();	}
		else	{
			$two_only++;	}
		$two=<FILE2>;
		($id_2, $flag_2, $chr_from_2, $loc_from_2, $mapq_2, $cigar_2, $d1_2, $d2_2, $d3_2, $read_2, $read_qual_2, @rest_2) = split(/\t/,$two);
	}
	while(<FILE2>)	{
	$two_only++;	}	
}	

if(eof(FILE2))	{
	while(!eof(FILE1) && ($id_1 lt $id_2))	{
		if ($id_1 eq $id_2)	{
			combine();	}
		else	{
			$one_only++;	}
		$one=<FILE1>;
		($id_1, $flag_1, $chr_from_1, $loc_from_1, $mapq_1, $cigar_1, $d1_1, $d2_1, $d3_1, $read_1, $read_qual_1, @rest_1) = split(/\t/,$one);
	}
	while(<FILE1>)	{
	$one_only++;	}	
}

if ($id_1 eq $id_2)	{
	combine();	
}

print STDERR "### Read pairs with both reads present:\t\t$both.\n";
print STDERR "### Read pairs with only read one present:\t$one_only.\n";
print STDERR "### Read pairs with only read two present:\t$two_only.\n";

#SUBROUTINES=============================================================================================================================================
			
sub dec2bin {
	my $str = unpack("B32", pack("N", shift));
	return $str;
}

sub bin2dec {
	return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}

#function to comapre read name as samtool sort by name
sub strnum_cmp{
	my $i = 0;
	my $j = 0;
	my $name_a = $_[0];
	my $name_b = $_[1];
	while ($i < length($name_a) && $j < length($name_b)){
		if (looks_like_number(substr($name_a,$i,1)) && looks_like_number(substr($name_b,$j,1))){
			while (substr($name_a,$i,1) eq '0'){$i+=1;}
			while (substr($name_b,$j,1) eq '0'){$j+=1;}
			while (looks_like_number(substr($name_a,$i,1)) && looks_like_number(substr($name_b,$j,1)) && substr($name_a,$i,1) eq substr($name_b,$j,1)){
				$i+=1;
				$j+=1;
			}
			if (looks_like_number(substr($name_a,$i,1)) && looks_like_number(substr($name_b,$j,1))){
				my $k = 0;
				while (looks_like_number(substr($name_a,$i+$k,1)) && looks_like_number(substr($name_b,$j+$k,1))) {$k+=1;}
				return ($i+$k < length($name_a) && looks_like_number(substr($name_a,$i+$k,1))) ? 1 : ($j+$k < length($name_b) && looks_like_number(substr($name_b,$j+$k,1))) ? -1 : int(substr($name_a,$i,1)) - int(substr($name_b,$j,1));
			}
			elsif (looks_like_number(substr($name_a,$i,1))) {return 1;}
			elsif (looks_like_number(substr($name_b,$j,1))) {return -1;}
			elsif ($i != $j) {return $i < $j ? 1 : -1;}
		}
		else{
			if (substr($name_a,$i,1) ne substr($name_b,$j,1)){
				return int(substr($name_a,$i,1)) - int(substr($name_b,$j,1));
			}
			$i+=1;
			$j+=1;
		}
	}
	return substr($name_a,$i,1) ? 1 : substr($name_b,$j,1) ?  -1 : 0;
}

sub lessthan{
	my $t = strnum_cmp($_[0],$_[1]);
	if ($t < 0){
		return 1;
	}
	else{
		return 0;
	}
}

sub combine	{
	my $bin1 = reverse(dec2bin($flag_1));
	my $bin2 = reverse(dec2bin($flag_2));

	my @binary1 = split(//,$bin1);
	my @binary2 = split(//,$bin2);

	if (($binary1[2] == 1) && ($mapq_1 >= 10)) {
		$mapq_1 = 0;	}
	if (($binary2[2]== 1) && ($mapq_2 >= 10)) {
		$mapq_2 = 0;	}

	my $proper_pair1;
	my $proper_pair2;
	my $dist1;
	my $dist2;

	if (($binary1[2] == 0) && ($binary2[2] == 0)) {
		$proper_pair1 = 1;
		$proper_pair2 = 1;
		if ($chr_from_1 eq $chr_from_2) {
			my $dist = abs($loc_from_1 - $loc_from_2);
			if ($loc_from_1 >= $loc_from_2) {
				$dist1 = -1*$dist;
				$dist2 = $dist;
			} else {
				$dist1 = $dist;
				$dist2 = -1*$dist;
			}
		} else {
			$dist1 = 0;
			$dist2 = 0;
		}
	} 
	else {
		$proper_pair1 = 0;	
		$proper_pair2 = 0;
		$dist1 = 0;
		$dist2 = 0;
	}

	my $new_bin1 = join("","000000000000000000000",$binary1[10],$binary1[9],$binary1[8],"0","1",$binary2[4],$binary1[4],$binary2[2],$binary1[2],$proper_pair1,"1");
	my $new_bin2 = join("","000000000000000000000",$binary2[10],$binary2[9],$binary2[8],"1","0",$binary1[4],$binary2[4],$binary1[2],$binary2[2],$proper_pair2,"1");

	my $new_flag1 = bin2dec($new_bin1);
	my $new_flag2 = bin2dec($new_bin2);

	print(join("\t",$id_1,$new_flag1,$chr_from_1,$loc_from_1,$mapq_1,$cigar_1,$chr_from_2,$loc_from_2,$dist1,$read_1,$read_qual_1,@rest_1));
	print(join("\t",$id_2,$new_flag2,$chr_from_2,$loc_from_2,$mapq_2,$cigar_2,$chr_from_1,$loc_from_1,$dist2,$read_2,$read_qual_2,@rest_2));
	$both++;
}

=BEGIN
my $bin1 = reverse(dec2bin($flag_1));
my $bin2 = reverse(dec2bin($flag_2));

my @binary1 = split(//,$bin1);
my @binary2 = split(//,$bin2);

if (($binary1[2] == 1) && ($mapq_1 >= 10)) {
	$mapq_1 = 0;	}
if (($binary2[2]== 1) && ($mapq_2 >= 10)) {
	$mapq_2 = 0;	}

my $proper_pair1;
my $proper_pair2;
my $dist1;
my $dist2;

if (($binary1[2] == 0) && ($binary2[2] == 0)) {
	$proper_pair1 = 1;
	$proper_pair2 = 1;
	if ($chr_from_1 eq $chr_from_2) {
		my $dist = abs($loc_from_1 - $loc_from_2);
		if ($loc_from_1 >= $loc_from_2) {
			$dist1 = -1*$dist;
			$dist2 = $dist;
		} else {
			$dist1 = $dist;
			$dist2 = -1*$dist;
		}
	} else {
		$dist1 = 0;
		$dist2 = 0;
	}
} 
else {
	$proper_pair1 = 0;
	$proper_pair2 = 0;
	$dist1 = 0;
	$dist2 = 0;
}

my $new_bin1 = join("","000000000000000000000",$binary1[10],$binary1[9],$binary1[8],"0","1",$binary2[4],$binary1[4],$binary2[2],$binary1[2],$proper_pair1,"1");
my $new_bin2 = join("","000000000000000000000",$binary2[10],$binary2[9],$binary2[8],"1","0",$binary1[4],$binary2[4],$binary1[2],$binary2[2],$proper_pair2,"1");

my $new_flag1 = bin2dec($new_bin1);
my $new_flag2 = bin2dec($new_bin2);

print(join("\t",$id_1,$new_flag1,$chr_from_1,$loc_from_1,$mapq_1,$cigar_1,$chr_from_2,$loc_from_2,$dist1,$read_1,$read_qual_1,@rest_1));
print(join("\t",$id_2,$new_flag2,$chr_from_2,$loc_from_2,$mapq_2,$cigar_2,$chr_from_1,$loc_from_1,$dist2,$read_2,$read_qual_2,@rest_2));
$both++;	

#if(eof(FILE2))	{
#	last;	}
=END
