##script to create expression profiles from a single sample

$usage = "perl expressionProfile_singleSample.pl input.patman output";
$inp = shift or die $usage;
$out = shift or die $usage;

##print "input file: $inp\n";
##print "output file: $out\n";

open inp, $inp
	or die "cannot open $inp";
open out, ">".$out
	or die "cannot open $out";

##initialize the gene name and the gene length
$geneName   = "";
$geneLength = 0;
for($j = 0; $j <= 50000; $j++)
{$expression[$j] = 0}

while(<inp>)
{
	chomp; 
	@d  = split(/\t/);
	@d0 = split(/ /,$d[0]);
	$d[0] = $d0[0];
	@d1 = split(/\-/,$d[1]);

#	print "#$d[0]# eq #$geneName#\n";
#	sleep(1);
	
	if($d[0] eq $geneName)
	{
		if($geneLength < $d[3]) {$geneLength = $d[3]}
		for($j = $d[2]; $j <= $d[3]; $j++)
		{
			$expression[$j] = $expression[$j]+$d1[1];
		}
	}
	else
	{
		if($geneLength > 0)
		{
			print out "$geneName, $geneLength,";
			$myString = "";
			for($j = 0; $j <= $geneLength; $j++)
			{
				$myString = $myString." ".$expression[$j];
			}
			print out "$myString\n";
		}
			$geneName = $d[0];
			$geneLength = 0;
			for($j = 0; $j <= 50000; $j++)
			{$expression[$j] = 0}
			
			if($geneLength < $d[3]) {$geneLength = $d[3]}
			for($j = $d[2]; $j <= $d[3]; $j++)
			{
				$expression[$j] = $expression[$j]+$d1[1];
			}
	}
}
close(inp);
close(out);

exit;
