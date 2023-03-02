
# This script has been written by Francesco Comandatore (Univerisyt of Milan)

# The scope of the the script is ti parallelize the script Infer_ancestral_state_changes.R

# For help run

# perl PARALLELIZE_Infer_ancestral_state_changes.pl

use POSIX;

# INPUTs

$pa1=shift;
$sh1=shift;
$pa2=shift;
$sh2=shift;
$pa3=shift;
$sh3=shift;
$pa4=shift;
$sh4=shift;

$hash{$pa1}=$sh1;
$hash{$pa2}=$sh2;
$hash{$pa3}=$sh3;
$hash{$pa4}=$sh4;

my $tab = $hash{'-tab'}; 	
my $tree=$hash{'-tree'};
my $num = $hash{'-cpu'};
my $prob = $hash{'-p'};

$Rscript_path = "Infer_ancestral_state_changes.R";

if ($prob eq ""){$prob = 0.75;}
if ($num eq ""){$num = 2;}

if ($tab eq "" or $tree eq "")
{
	print "This script allows to infer the ancestral state of each character present in table in each node of a phylogenetic tree. Bootstrap supports are not considered.

Run: perl script.pl -tab [table] -tree [tree] -p [probability for ancestral reconstruction] -cpu [cpu]

TABLE:

Example of table (TAB delimited):

ID_SNP	feature1	feature2	org1	org2 org3
SNP_1	Synonimous	Gene annotation	a	a	t
SNP_2	Synonimous	Gene annotation	c	g	g
 

notes: ID_SNP must be unique, does not matter the number or position of feature columns

TREE:

((org1,org2),org3)

All the organism present in the tree must be present in the table!

";

exit;

}


##############################################################################################################à

open(TTT,"<$tab");
$title=<TTT>;
close(TTT);

`sed 1d $tab > $tab.nohead`;
`sort -k4 $tab.nohead > $tab\_sorted`;
`rm $tab.nohead`;

open(TTT,"<$tab\_sorted");

$num_row=`wc -l $tab\_sorted`;
$div=ceil($num_row/$num)+1;

mkdir("tmp");

$c=1;
$file=1;
open(OUT,">tmp/$tab\_sorted.$file");
print OUT $title;

while(<TTT>)
{
	print OUT $_;

	$c++;

	if ($c == $div)
	{close(OUT);$file++;open(OUT,">tmp/$tab\_sorted.$file");print OUT $title;$c=1;}	

}

close(OUT);

if ($num_row % $num == 0){`rm tmp/$tab\_sorted.$file`;}

# Run Rscript

`ls tmp/* | parallel "Rscript $Rscript_path {} $tree $prob"`;

# Get header

@ls_tmp = `ls tmp/*.anc_state_changes.tab`;
chomp $ls_tmp[0];

`head -n1 $ls_tmp[0] > $tab.anc_state_changes.tab`;

# Get all the lines from tmp files

`sed -i 1d tmp/*.anc_state_changes.tab`;

`cat tmp/*.anc_state_changes.tab >> $tab.anc_state_changes.tab`;

`rm $tab\_sorted`;
`rm -r tmp`;


