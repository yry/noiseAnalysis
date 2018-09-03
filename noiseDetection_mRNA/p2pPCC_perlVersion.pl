##input csv file with expression levels
##ouput p2p PCC and abundance matrices
##in the output matrix on position [i,k] we represent the information for sample i, transcript k

$usage = "perl p2pPCC_perl.pl expressionLevel_matrix.csv sampleCount PCCmatrix ABNmatrix";
$inp   = shift or die $usage; ##input in the pre-defined format
$sCt   = shift or die $usage; ##number of samples to be analysed
$pcc   = shift or die $usage; ##matrix of p2p PCCs
$abn   = shift or die $usage; ##matrix of abundances

open in, $inp
	or die "cannot open input";
open pcc, ">".$pcc
	or die "cannot open pcc output";
open abn, ">".$abn
	or die "cannot open abn output";
$transcript = "";
$transcriptLength = 0;

##initialise the profiles
for($i=0; $i <= $sCt; $i++)
{
	for($j = 0; $j < 100; $j++)
	{$profiles[$i][$j] = 0;}
}

while(<in>)
{
	chomp;
	@d =split(/,/);
	## first retrieve the profiles corresponding to a transcript
	if($transcript eq $d[1])
	{
		@dd = split(/ /,$d[3]);
		#$d[2] is the transcript length
		for($k=0; $k < $d[2]; $k++)
		{$profiles[$d[0]][$k] = $dd[$k];}
		$transcriptLength = $d[2];
	}
	else
	{
		if($transcript ne "")
		{
			#$abnVector will be the vector of abundances 
			#$pccVector will be the vector of p2pPCC values
			#print "working with transcript $transcript of length $transcriptLength\n";
			##calculate the PCC
			##calculate the maxY abundance
			for($j = 0; $j < $sCt; $j++)
			{
				$maxY = 0;
				for($k = 0; $k < $transcriptLength; $k++)
				{if($maxY < $profiles[$j][$k]){$maxY = $profiles[$j][$k];}}
				#print "maxY is $maxY\n";
				
				if($maxY > 1)
				{
				##calculate the PCCs
				$meanSample = 0;
				for($k = 0; $k < $transcriptLength; $k++)
				{
					$meanSample = $meanSample + $profiles[$j][$k];
				}
				$meanSample = $meanSample/$transcriptLength;

				$#pcc= -1;				
				for($i = 0; $i < $sCt; $i++)
				{
				if($i != $j)
					{	
						$mean2 = 0;
						for($k = 0; $k < $transcriptLength; $k++)
						{
							$mean2 = $mean2 + $profiles[$i][$k];
						}
						$mean2 = $mean2/$transcriptLength;

						$meanSample  = int($meanSample * 100) / 100;
						$mean2 = int($mean2 * 100) / 100;
						
						#print "mean sample: $meanSample\n";
						#print "for sample $i the mean is $mean2\n";
						
						if($mean2 == 0)
						{
							$r = 0;
						}
						else
						{
							$numerator = 0;$denominator1=0;$denominator2=0;
							for($k = 0; $k < $transcriptLength; $k++)
							{
								$numerator = $numerator + ($profiles[$j][$k]-$meanSample)*($profiles[$i][$k]-$mean2);
								$denominator1 = $denominator1 + ($profiles[$j][$k]-$meanSample)*($profiles[$j][$k]-$meanSample);
								$denominator2 = $denominator2 + ($profiles[$i][$k]-$mean2)*($profiles[$i][$k]-$mean2);
							}
							$numerator = int($numerator * 100) / 100;
							$denominator1 = int($denominator1 * 100) / 100;
							$denominator2 = int($denominator2 * 100) / 100;
							
							$r = $numerator/(sqrt($denominator1)*sqrt($denominator2));
							$r = int($r * 100) / 100;
						}
						#print "num = $numerator, den1: $denominator1, den2: $denominator2\n";
						#sleep(1);
						#print "with sample $i the pcc is $r\n";
						
						push @pcc, [($r)];
					
					}##endif
				}##end for i
				$meanPCC = 0;
				for($i=0; $i <= $#pcc; $i++)
				{
					#print "pcc[$i] is $pcc[$i][0]\n";
					$meanPCC = $meanPCC + $pcc[$i][0];
				}
				$meanPCC = $meanPCC / ($#pcc + 1);

				$abnVector[$j] = $maxY;
				$pccVector[$j] = $meanPCC;
				#print "sample $j with abn $maxY and average PCC $meanPCC\n";
				#sleep(1);

				}##end for maxY > 1
				else
				{
					$abnVector[$j] = 0;
					$pccVector[$j] = 0;
				}
			}
			for($j = 0; $j <$sCt-1; $j++)
			{
				print pcc "$pccVector[$j],";
				print abn "$abnVector[$j],";
			}
			print pcc "$pccVector[$j]\n";
			print abn "$abnVector[$j]\n";
		}##end of the PCC calculations and prints
		
		$transcript = $d[1];
		##reinitialize the profiles
		for($i=0; $i <= $sCt; $i++)
		{
			for($j = 0; $j < 100; $j++)
			{$profiles[$i][$j] = 0;}
		}
		@dd = split(/ /,$d[3]);
		#$d[2] is the transcript length
		for($k=0; $k < $d[2]; $k++)
		{$profiles[$d[0]][$k] = $dd[$k];}
		$transcriptLength = $d[2];
	}

}
close(in);
close(pcc);
close(abn);

exit;

