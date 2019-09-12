#!/usr/bin/perl

use Cwd 'abs_path';
use strict;
use warnings;
use Getopt::Long qw(:config pass_through);

my ($infile, $taxa, $outfile, $filter) = ("", "", "", "F");
GetOptions (
  "i|group:s" => \$infile,
  "t|taxon:s" => \$taxa,
  "o|out:s" => \$outfile,
  "f|filter:s" => \$filter,
);

print("Open file $taxa\n");
open(F1,$taxa) or die "Cant open : $! \n";

my $OUTFILE;
open $OUTFILE, '>', $outfile
    	or die "Cannot open $outfile";

#use a hash hash to save the data with the structure:
#
#  first_hash   /---->key2=parent_id,value=parent taxid
# key1=taxid---{
#               \---->key3=tax_name, value=taxon name
#                          second_hash
#
################################################################
# open taxonomy file and save it in hashname with the hash structure mentioned above

my %hashname;

my $pa="parent_id";
my $ta="tax_name";

while (my $line =<F1> )  {
	chomp($line);
	$line =~ s/\"//g;
	my @arr = split(",", $line);
	$hashname{$arr[0]}{$pa}=$arr[1];
	$hashname{$arr[0]}{$ta}=$arr[3];
}

close(F1);

###################################################################

#open the ortholog group file

print("Open file $infile\n");
open(F2,$infile) or die "Cant open : $! \n";

my $userinput;
if($filter eq "T"){
print "Group:speciesA 1 sequence, speciesB __ sequences.\nGive the number of sequences of speciesB:\n";
$userinput =  <STDIN>;
chomp ($userinput);
}

my %parent;

while (my $line =<F2> )  {
	chomp($line);
	my @arr = split(" ", $line);
	#print("$arr[0]\n$arr[1]\n$arr[2]\n");
	print{$OUTFILE}("$arr[0]\t");
	my %tax;
	for(my $i=1;$i<=$#arr;$i++){
		$arr[$i] =~ /(\d+)\|.+/g;
		
		if(exists($tax{$1})){
			$tax{$1}+=1;
			#print("$1\t$tax{$1}\n");
		}else{
			$tax{$1}=1;
		}
	}

# if the group is like: speciesA 1sequence, speciesB 1000 sequences, then the common ancestor is speciesB, consider speciesA as outlier.
# user can set the threshold, speciesB can be 1000 sequeces, or 500 sequeces....

	if($filter eq "T"){
#		print "Group:speciesA 1 sequence, speciesB __ sequences.\nGive the number of sequences of speciesB:\n";
#		my $userinput =  <STDIN>;
#		chomp ($userinput);
		my $ttt=0;
		my $tta="";
		my $ttb=0;
		if(keys %tax ==2){
			for (keys %tax){
				if($tax{$_}==1){
					$ttt +=1;
					$tta = $_;
				}
				if($tax{$_}>=$userinput){
					$ttb = 1;
				}
			}
			if(($ttt==1)&&($ttb==1)){
				delete $tax{$tta};	
			}
		}
	}

	my @ke = keys %tax;
	#foreach(@ke){
#print("$_\n");}
	if(@ke==1){
		print{$OUTFILE}("$ke[0]\t$hashname{$ke[0]}{$ta}\t1\t$ke[0]\t$tax{$ke[0]}\n");
	}else{
		my $co = CommonAncestor(@ke);
		#print("Find common ancestor $co \n");
		my $num=$#ke+1;
		print{$OUTFILE}("$co\t$hashname{$co}{$ta}\t$num");
############################################################################
# To print the ancestor you are interested:
# print{$OUTFILE}("\t$hashname{4751}{$ta}\t$parent{4751}");
# change 4751 to the taxid you like.
############################################################################
		if($co==33213){print{$OUTFILE}("\t$hashname{33511}{$ta}\t$parent{33511}\t$hashname{1206794}{$ta}\t$parent{1206794}\t$hashname{1206795}{$ta}\t$parent{1206795}");}
		if($co==33154){print{$OUTFILE}("\t$hashname{4751}{$ta}\t$parent{4751}\t$hashname{33208}{$ta}\t$parent{33208}\t$hashname{10226}{$ta}\t$parent{10226}\t$hashname{6040}{$ta}\t$parent{6040}\t$hashname{33213}{$ta}\t$parent{33213}\t$hashname{33511}{$ta}\t$parent{33511}\t$hashname{1206794}{$ta}\t$parent{1206794}\t$hashname{1206795}{$ta}\t$parent{1206795}\t$hashname{6073}{$ta}\t$parent{6073}\t$hashname{10197}{$ta}\t$parent{10197}");}
		if($co==6072){print{$OUTFILE}("\t$hashname{33213}{$ta}\t$parent{33213}\t$hashname{33511}{$ta}\t$parent{33511}\t$hashname{1206794}{$ta}\t$parent{1206794}\t$hashname{1206795}{$ta}\t$parent{1206795}\t$hashname{6073}{$ta}\t$parent{6073}\t$hashname{10197}{$ta}\t$parent{10197}");}
		if($co==33208){print{$OUTFILE}("\t$hashname{10226}{$ta}\t$parent{10226}\t$hashname{6040}{$ta}\t$parent{6040}\t$hashname{33213}{$ta}\t$parent{33213}\t$hashname{33511}{$ta}\t$parent{33511}\t$hashname{1206794}{$ta}\t$parent{1206794}\t$hashname{1206795}{$ta}\t$parent{1206795}\t$hashname{6073}{$ta}\t$parent{6073}\t$hashname{10197}{$ta}\t$parent{10197}");}
		my $sum;

		foreach my $k (keys %tax){
			$sum+=$tax{$k};
			print{$OUTFILE}("\t$k\t$tax{$k}");
		}
		print{$OUTFILE}("\tSum: $sum\n");
	}
	%parent=();
}

sub CommonAncestor{
	my @list=@_;
	my $l=$#list+1;
	my @children=@_;
	#my %parent;
	my %rank;
	my $com;
	my $m = 1;
	my $r = 0;
	do{
	my @par;
	foreach my $id (@list){
		if(exists($parent{$id})){
			$parent{$id}+=1;
			$rank{$id} = $r;
	#	print("$id\t$parent{$id}\n");
		}else{
			$parent{$id} = $m;
			$rank{$id} = $r;
		}
		push(@par,$hashname{$id}{$pa});
		}
	@list=@par;
	$r+=1;
	}while((!exists($parent{1}))||($parent{1}<100*$l));

        my @keys = sort { $rank{$a} <=> $rank{$b}} keys(%rank);
        for(my $k=0;$k<=$#keys;$k++){
        #print("$keys[$k]\t$parent{$keys[$k]}\n");              
        if($parent{$keys[$k]}==$l){
                        $com=$keys[$k];
                        return($com);
                }
        }

	#return($com);
}
