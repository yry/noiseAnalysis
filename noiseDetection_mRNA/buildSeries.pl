##group the expression profiles from the different samples
##collate all single sample expression profiles
##this will be intput for the p2p PCC

sub trim($);
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

$usage = "perl build series.pl out";
$out = shift or die $usage;

opendir(DIR, ".");
@files = grep(/\_final.expression$/,readdir(DIR));
closedir(DIR);

print "number of samples to be combine: $#files\n";

##maintain a file with the order of the samples
open cc, ">combinedOrder.csv"
	or die "cannot open combinedOrder.csv";

$i = 0;
foreach $file (@files) 
{
	print cc "$file\n";
	print "working with $file\n";
#	system("date");
	
	open f, $file
		or die "cannot open $file\n";
	while(<f>)
	{
		chomp;
		@d = split(/,/);
		push @profiles, [($i,trim($d[0]),trim($d[1]),trim($d[2]))];
	}
	close(f);
	print "profiles read in sample $i [$file]: $#profiles\n";
	$i++;
}
close(cc);

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
#	print "0: $ps[$j][0]\n";
#	print "1: $ps[$j][1]\n";
#	print "2: $ps[$j][2]\n";
#	print "3: $ps[$j][3]\n";
#	print "memo: $geneMemo\n";
#	sleep(1);
	
	if($ps[$j][1] eq $geneMemo)
	{
		push @combine, [($ps[$j][0],$ps[$j][1],$ps[$j][2],$ps[$j][3])];
		$completion[$ps[$j][0]] = 1;
		$completedSamples++;
#		print "completed samples: $completedSamples and entries in combined: $#combine\n";
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
					push @combine, [($k, $ps[$j-1][1], 1, "1 0 ")];
				}
			}
		}
		##sort on the sample identifier, to have the components in the same order
		@cs = sort{$a->[0] <=> $b->[0]}@combine;
		
		##process the expression profiles, to ensure the same length on all expression profiles
		$maxLen = 0;
		for($k = 0; $k <= $#files; $k++)
		{
			if($cs[$k][2] > $maxLen){$maxLen = $cs[$k][2]}
		}
		
#		print "overall max length is: $maxLen and final number of entries: $#cs\n";
		
		for($k = 0; $k <= $#files; $k++)
		{
			print out "$cs[$k][0], $cs[$k][1], $maxLen, $cs[$k][3]";
			if($cs[$k][2] < $maxLen)
			{
				##0 padding for shorter versions of the gene
				for($l = $cs[$k][2]; $l <= $maxLen; $l++)
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
		
		$#combine = -1;
		push @combine, [($ps[$j][0],$ps[$j][1],$ps[$j][2],$ps[$j][3])];
		$completion[$ps[$j][0]] = 1;
		$completedSamples++;
	}
	
}
close(out);

exit;