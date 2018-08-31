##get the file identity for the combined files
##invoke from a folder with csv files
##provide the corresponding fasta file

use Bio::SeqIO;

sub trim($);
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

$usage = "perl getFileIdentity.pl fasta.out out";
$faOut = shift or die $usage;
$out   = shift or die $usage;

##
## read and organize the csv files
##

opendir(DIR, ".");
@files = grep(/\.csv$/,readdir(DIR));
closedir(DIR);

$index = 1;
foreach $f (@files) 
{
	chomp;
	$fileName = $f;
	print "working with: $fileName\n";
	open in, $fileName
		or die "cannot open file $fileName";
	while(<in>)
	{
		chomp;
		@d = split/,/;
		##if($d[1] >= 2) ## this was the condition when the files were created.
		{
		##sequence, abundance, and index of the file
		push @gIn, [($d[0], $d[1],$index)];
		}
	}
	close(in);
	$index++;
	print "entries In: $#gIn\n";
}
close(inp);
$index--;

print "files read: $index and entries: $#gIn \n";

@gs = sort{$a->[0] cmp $b->[0]}@gIn;
print "sorting finished. \n";

$select = 1; 
#1 is for R expression
#2 is for NR expression

$gene = $gs[0][0];
for($j = 1; $j <= $index; $j++)
{
	$expr[$j] = 0;
}
$expr[$gs[0][3]] = $gs[0][$select];

for($i = 1; $i <= $#gs; $i++)
{
	if($gs[$i][0] eq $gene)
	{
		$expr[$gs[$i][2]] = $gs[$i][$select];
	}
	else
	{
		$expr[0] = $gene;
		$mySum = 0;
		for($j = 1; $j <= $index; $j++)
		{$mySum = $mySum + $expr[$j]}

#		print "mySum: $mySum and expr[0]: $expr[0]\n";
#		sleep (2);

		if($mySum > 2 && length($expr[0]) > 15)
		{
		push @expression, [@expr];}
		$gene = $gs[$i][0];
		for($j = 1; $j <= $index; $j++)
		{
			$expr[$j] = 0;
		}
		$expr[$gs[$i][2]] = $gs[$i][$select];
	}
}
$expr[0] = $gene;
$mySum = 0;
for($j = 1; $j <= $index; $j++)
{$mySum = $mySum + $expr[$j]}
if($mySum > 2 && length($expr[0]) > 15)
{push @expression, [@expr];}

print "expression series created: $#expression \n";

##
## read the input fasta files
##
my $fa = Bio::SeqIO->new('-format' => 'fasta', '-file' =>  "$faOut");
while(my $seq = $fa->next_seq()) 
{
	my $sequence = $seq->seq();
	my $id       = $seq->display_id();
	@idd = split(/_/,$id);
	push @faExpr, [(@idd)];
}
print "fasta file read, with $#faExpr entries.\n";

##
## compare the file from the csv files and the file from the fasta
##

print "number of samples from fasta: $#idd; number of rows $#faExpr\n";
print "number of samples from the csv: $index; number of rows $#expression\n";
$nrSamples = $index;
#sleep(100);

open out, ">".$out
	or die "cannot open output";

##the fasta matrix is the pivot comparator
for($j = 1; $j <= $nrSamples; $j++)
{
	print "working with column $j\n";
	$selectedSample = -1; ## flag variable whether a sample is the selected one
	for($k = 1; $k <= $nrSamples; $k++)
	{
		$okSample = 1;
		for($i = 0; $i <= $#faExpr; $i++)
		{
			if($faExpr[$i][$j] != $expression[$i][$k])
			{$okSample = -1; $i = $#faExpr;}
		}
		if($okSample == 1)
		{$selectedSample = $k; $k = $nrSamples + 1;}
	}
	print out "$files[$selectedSample-1]\n";
	print "$j and $files[$selectedSample-1]\n";
}
close(out);

exit;

