##script for processing large files with reads in non-redundant format

use Bio::SeqIO;

$usage = "perl createExpressionProfiles_largeNRfiles.pl transcriptomicsReference";
$genomeReference = shift or die $usage;

$batchSize = 1000000;

opendir(DIR, ".");
@files = grep(/\.fasta$/,readdir(DIR));
closedir(DIR);

foreach $file (@files) 
{
	print "working with $file\n";
	system("date");
	
	@f = split(/\./,$file);
	$patmFile = $f[0].".patman";
	$exprFile = $f[0].".expression";
	$exprFileFinal = $f[0]."_final.expression";
	print "patman file: $patmFile\n";
	print "expression file: $exprFile\n";
	
	my $fa = Bio::SeqIO->new('-format' => 'fasta', '-file' =>  "$file");
	while(my $seq = $fa->next_seq()) 
	{
		my $sequence = $seq->seq();        # get the sequence as a string
		my $id       = $seq->display_id();  # identifier
		@idd = split(/\_/,$id);
		$newID = $idd[0]."-".$idd[1];
		push @seqs, [$newID, $sequence];
	}
	print "number of entries in the file: $#seqs\n";
	$batches = int($#seqs / $batchSize);
	
	print "the reads will be processed in $batches batches\n";
	
	$exprFiles = ""; ##create a string with the list of expression files
	for($k = 0; $k <= $batches; $k++)
	{
		open outTmp, ">outTMP.fasta"
			or die "cannot open outTMP on batch $k";
		$exprFileK = $exprFile."_$k";
		for($j = 0; $j < $batchSize; $j++)
		{
			if(length($seqs[$k*$batchSize+$j][1]) > 0)
			{
				print outTmp ">$seqs[$k*$batchSize+$j][0]\n$seqs[$k*$batchSize+$j][1]\n";
			}
		}
		#print "patman -e 0 -g 0 -D $genomeReference -P outTMP.fasta -o patmTMP.patman\n";
		system("patman -e 0 -g 0 -D $genomeReference -P outTMP.fasta -o patmTMP.patman");
		print "perl expressionProfile_singleSample.pl patmTMP.patman $exprFileK\n";
		system("perl expressionProfile_singleSample.pl patmTMP.patman $exprFileK");
		system("rm patmTMP.patman");
		$exprFiles = $exprFileK." ".$exprFiles;
	}
	
	##merge the expression files
	print "cat $exprFiles > $exprFile\n";
	system("cat $exprFiles > $exprFile");
	
	##remove the intermediary files
	for($k = 0; $k <= $batches; $k++)
	{
		$exprFileK = $exprFile."_$k";
		system("rm $exprFileK");
	}
	
	##process the expression file such that all entries/genes apprea only once
	print "expression file: $exprFile\n";
	open expr, $exprFile
		or die "cannot open expression file";
	open out, ">".$exprFileFinal
		or die "cannot open $exprFileFinal";
		
	while(<expr>)
	{
		chomp;
		@d = split(/,/);
		push @expressionProfiles, [@d];
	}
	close(expr);
	print "expression profiles read: $#expressionProfiles\n";
	
	@es = sort{uc($a->[0]) cmp uc($b->[0])}@expressionProfiles;
	print "sorting completed.\n";
	
	$geneMemo = $es[0][0];
	for($i = 0; $i <= $#es; $i++)
	{
		if($geneMemo eq $es[$i][0])
		{
			push @combine,[($es[$i][0],$es[$i][1],$es[$i][2])];
		}
		else
		{
			$maxLen = $combine[0][1];
			for($k = 0; $k <= $#combine; $k++)
			{
				if($combine[$k][1] > $maxLen)
				{$maxLen = $combine[$k][1]}
			}
			##create the combined expression levels
			for($k = 0; $k <= $#combine; $k++)
			{
				@c = split(/ /,$combine[$k][2]);
				for($l = 0; $l <= $maxLen; $l++)
				{
					$expression[$l] = $expression[$l] + $c[$l];
				}
			}
			
			$myString = $expression[0];
			for($l = 0; $l <= $maxLen; $l++)
			{
				$myString = $myString." ".$expression[$l];
			}
			print out "$geneMemo, $maxLen, $myString\n";
			
			##reinitialize
			for($l = 0; $l <= 10000; $l++)
			{$expression[$l] = 0}
			
			##restart with the new entry
			$geneMemo = $es[$i][0];
			$#combine = -1;
			push @combine,[($es[$i][0],$es[$i][1],$es[$i][2])];
		}
	}
	close(out);
} ##end of foreach


exit;