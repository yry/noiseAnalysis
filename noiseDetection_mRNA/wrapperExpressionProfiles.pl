##wrapper expression profiles

$usage = "wrapperExpressionProfiles.pl transcripts";
$transcripts = shift or die $usage;

system("perl createExpressionProfiles_largeNRfiles.pl $transcripts");
system("perl buildSeries.pl allSeries.expression");

opendir(DIR, ".");
@files = grep(/\final.expression$/,readdir(DIR));
closedir(DIR);
print "number of samples to work with: $#files\n";

system("perl p2pPCC_perlVersion.pl allSeries.expression $#files matrix.pcc matrix.abn");

exit;