#combine expression levels for genes

sub trim($);
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

$usage = "perl combine.pl out";
$out = shift or die $usage;

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
		push @gIn, [($d[0], $d[1], $d[2],$index)];
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
		$expr[$gs[$i][3]] = $gs[$i][$select];
	}
	else
	{
		$expr[0] = $gene;
		push @expression, [@expr];
		$gene = $gs[$i][0];
		for($j = 1; $j <= $index; $j++)
		{
			$expr[$j] = 0;
		}
		$expr[$gs[$i][3]] = $gs[$i][$select];
	}
}
$expr[0] = $gene;
push @expression, [@expr];

print "expression series created: $#expr \n";

open outCSV, ">".$out.".csv"
	or die "cannot open output";
open outFASTA, ">".$out.".fasta"
	or die "cannot open output";

for($i = 0; $i <= $#expression; $i++)
{
	$id = "";$mySum = 0;
	for($j = 1; $j <= $index; $j++)
	{$mySum = $mySum + $expression[$i][$j]}

	if($mySum > 2 && length($expression[$i][0]) > 15)
	{
		for($j = 0; $j < $index; $j++)
		{
			$expression[$i][$j] = trim($expression[$i][$j]);
			print outCSV "$expression[$i][$j], ";
			$id = $id.$expression[$i][$j]."_";
		}
		$expression[$i][$j] = trim($expression[$i][$j]);
		print outCSV "$expression[$i][$j]\n";
		$id = $id.$expression[$i][$j];

		print outFASTA ">$id\n$expression[$i][0]\n";
	}
}
close(outCSV);
close(outFASTA);
exit;
