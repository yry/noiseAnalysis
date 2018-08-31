##group the expression profiles from the different samples
##collate all single sample expression profiles
##this will be intput for the p2p PCC

$usage = "perl build series.pl out";
$out = shift or die $usage;

opendir(DIR, ".");
@files = grep(/\final.expression$/,readdir(DIR));
closedir(DIR);

print "number of samples to be combine: $#files\n";

$i = 0;
foreach $file (@files) 
{
	print "working with $file\n";
	system("date");
	
	open f, $file
		or die "cannot open $file\n";
	while(<f>)
	{
		chomp;
		@d = split(/,/);
		push @profiles, [($i,$d[0],$d[1],$d[3])];
	}
	close(f);
	print "profiles read in sample $i [$file]: $#profiles\n";
	$i++;
}

##sort on the gene identifier (string value) and on the sample identifier
@ps = sort{uc($a->[1]) cmp uc($b->[1]) || $a->[0] <=> $b->[0]}@profiles;
print "sorting done.\n";

open out, ">".$out
	or die "cannot open $out";

##set of flags to ensure that for every gene we have the full set of expression levels
##if a gene is not present in a sample, that a full 0 expression will be introduced
##put a 1 value on the first position, to avoi 0 variance on the calculation of the Pearson Correlation Coeeficient

$geneMemo = $ps[0][1];
for($k = 0; $k < $#files; $k++)
{$completion[$k] = 0;}
$completedSamples = 0;

for($j = 0; $j <= $#ps; $j++)
{
	if($ps[$j][1] eq $geneMemo)
	{
		push @combine, [($ps[$j][0],$ps[$j][1],$ps[$j][2],$ps[$j][3])];
		$completion[$ps[$j][0]] = 1;
		$completedSamples++;
	}
	else
	{
		if($completedSamples < $#files)
		{
			##find the missing samples and append them
			for($k = 0; $k <= $#files; $k++)
			{
				if($completion[$k] == 0)
				{
					## the [1] is the gene identifier
					push @combine, [($k, $ps[$j-1][1], 1, "1 0")];
				}
			}
			##sort on the sample identifier, to have the components in the same order
			@cs = sort{$a->[0] <=> $b->[0]}@combine;
		}
		
		##process the expression profiles, to ensure the same length on all expression profiles
		$maxLen = 0;
		for($k = 0; $k <= $#files; $k++)
		{
			if($combine[$k][2] > $maxLen){$maxLen = $combine[$k][2]}
		}
		for($k = 0; $k <= $#files; $k++)
		{
			print out "$combine[$k][0], $combine[$k][1], $maxLen, $combine[$k][3]";
			if($combine[$k][2] < $maxLen)
			{
				##0 padding for shorter versions of the gene
				for($l = $combine[$k][2]; $l <= $maxLen; $l++)
				{
				print out "0 "
				}
			}
			print out "\n";
		}
		
		##reinitialize
		$geneMemo = $ps[$j][1];
		for($k = 0; $k < $#files; $k++)
		{$completion[$k] = 0;}
		$completedSamples = 0;
		
		push @combine, [($ps[$j][0],$ps[$j][1],$ps[$j][2],$ps[$j][3])];
		$completion[$ps[$j][0]] = 1;
		$completedSamples++;
	}
	
}
close(out);

exit;