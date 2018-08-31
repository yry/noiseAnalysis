##script to determine the expression profiles, per sample
##these outputs will then be comcatenated (all samples combine)
##this version of the script is for handling concatenated sets
##and the p2pPCC will be computed between samples

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

$usage = "perl noise_analysis_per_sample.pl genome winLen fasta_nr nrSamples mm out";
$genome = shift or die $usage; ##genome file
$win    = shift or die $usage; ##window length for the genome processing
$fastaNR= shift or die $usage; ##input file for which the expression profiles are calculated
$nrSamples = shift or die $usage; ##number of samples in fasta nr file
$mm     = shift or die $usage; ##number of mis-matches
$out    = shift or die $usage;

#@f      = split(/\./,$fastaNR);
#$out    = $f[0];
#for($k = 1; $k < $#f; $k++)
#{
#	$out = $out."_".$f[$k];
#}
#$out   = $out.".expression";
print "output file: $out\n";

my $fa = Bio::SeqIO->new('-format' => 'fasta', '-file' =>  "$genome");
while(my $seq = $fa->next_seq()) 
{
	my $sequence = $seq->seq();
	my $id       = $seq->display_id();
	push @gen, [($id,$sequence)];
}
print "genome read.\n";

$batch = 100;
$batchNo = int($#gen/$batch);

open out, ">".$out
	or die "cannot open output";
	
for($j = 0; $j < $batchNo; $j++)
{
	print "working with batch: $j\n";
	open oo, ">out_temp.fasta"
		or die "cannot open out_temp.fasta";
	for($i = $j * $batch; $i < ($j+1)*$batch; $i++)
	{
		print oo ">$gen[$i][0]\n$gen[$i][1]\n";
	}
	close(oo);
	system("patman -e $mm -g 0 -D out_temp.fasta -P $fastaNR -o gOut.patman");
	print "patman mapping finished\n";
	##process the temporary patman file
	
	$transcript = "";
	for($jj = 0; $jj <= $win; $jj++)
	{
		for($kk = 0; $kk < $nrSamples; $kk++)
		{$expr[$kk][$jj]=0}
	}

	open inp, "gOut.patman"
		or die "cannot open gOut.patman";
	while(<inp>)
	{
		chomp;
		@d = split(/\t/);
		for($kk = 0; $kk <= $#d; $kk++)
		{$d[$kk] = trim($d[$kk])}
		
		@dd= split(/\_/,$d[1]);
		
		if($d[0] eq $transcript)
		{
			for($kk=1; $kk <= $nrSamples; $kk++)
			{
				for($jj = $d[2]; $jj <= $d[3]; $jj++)
				{$expr[$kk-1][$jj] = $expr[$kk-1][$jj] + $dd[$kk]}
			}
		}
		else
		{
			if($transcript ne "")
			{
				for($kk=0; $kk < $nrSamples; $kk++)
				{
					print out "$kk,$transcript,$win,";
					$myString = "";
					for($jj = 0; $jj < $win-1; $jj++)
						{$myString = $myString."$expr[$kk][$jj] "}
					$myString = $myString."$expr[$kk][$jj]\n";
					print out "$myString";
				}
			}
			$transcript = $d[0];
			for($jj = 0; $jj <= $win; $jj++)
			{
				for($kk = 0; $kk < $nrSamples; $kk++)
				{$expr[$kk][$jj]=0}
			}
			
			for($kk=1; $kk <= $nrSamples; $kk++)
			{
			for($jj = $d[2]; $jj <= $d[3]; $jj++)
			{$expr[$kk-1][$jj] = $expr[$kk-1][$jj] + $dd[$kk]}
			}
		}
	}
	close(inp);
}
close(out);

##final clean up
system("rm out_temp.fasta");
system("rm gOut.patman");

exit;
