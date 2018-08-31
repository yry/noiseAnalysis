use Bio::SeqIO;

$usage = "perl noise_analysis_per_sample.pl genome winLen";
$genome = shift or die $usage; ##genome file
$win    = shift or die $usage; ##window length for the genome processing

my $fa = Bio::SeqIO->new('-format' => 'fasta', '-file' =>  "$genome");
while(my $seq = $fa->next_seq()) 
{
	my $sequence = $seq->seq();
	my $id       = $seq->display_id();
	if(index($id,"_") == -1)
	{
		print "id: $id\n";
		$out = $id."_win".$win.".fasta";
		$winCount = int(length($sequence)/$win)+1;
		open out, ">".$out
			or die "cannot open $out";
		print "number of windows to be created: $winCount\n";
		for($k = 0; $k < $winCount; $k++)
		{
			$seq = uc(substr($sequence,$k*$win, $win));
			$idf = $id."_".$k;
			if(index($seq,"N") == -1)
			{
			print out ">$idf\n$seq\n";}
		}
		close(out);
	}
}

exit;
