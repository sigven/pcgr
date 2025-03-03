use strict;


my $inputfile = $ARGV[0];
my $usemaxent = 1;


my @metables = &makemaxentscores;

open (FILE,"<$inputfile") || die "can't open!\n";

while(<FILE>) {
    chomp;
    if (/^\s*$/) { #discard blank lines;
	next;
    } 
    elsif (/^>/) { #discard comment lines;
	next;
    }
    elsif (/[NQWERYUIOPLKJHFDSZXVBM]/) {
	next;
    }
    else {
        $_ =~ s/\cM//g; #gets rid of carriage return
	my $str = $_;
	print $str."\t";
	$str = uc($str);
	if ($usemaxent) { 
	    print sprintf("%.2f", 
				     &log2(&scoreconsensus($str)*&maxentscore(&getrest($str),\@metables)))."\n";
	}
    }
}

sub hashseq{
    #returns hash of sequence in base 4
    # &hashseq('CAGAAGT') returns 4619
    my $seq = shift;
    $seq = uc($seq);
    $seq =~ tr/ACGT/0123/;
    my @seqa = split(//,$seq);
    my $sum = 0;
    my $len = length($seq);
    my @four = (1,4,16,64,256,1024,4096,16384);
    my $i=0;
    while ($i<$len) {
        $sum+= $seqa[$i] * $four[$len - $i -1] ;
	$i++;
    }
    return $sum;
}

sub makemaxentscores{
    my $dir = "splicemodels/";
    my @list = ('me2x3acc1','me2x3acc2','me2x3acc3','me2x3acc4',
		'me2x3acc5','me2x3acc6','me2x3acc7','me2x3acc8','me2x3acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	#print STDERR $file."\t".$num."\t".$n."\n";
	$num++;
    }
    return @metables;
}
sub makewmmscores{
    my $dir = "/bionet/geneyeo_essentials/MaxEntropy/webserver/splicemodels/";
    my @list = ('me1s0acc1','me1s0acc2','me1s0acc3','me1s0acc4',
		'me1s0acc5','me1s0acc6','me1s0acc7','me1s0acc8','me1s0acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	#print STDERR $file."\t".$num."\t".$n."\n";
	$num++;
    }
    return @metables;
}
sub makemmscores{
    my $dir = "/bionet/geneyeo_essentials/MaxEntropy/webserver/splicemodels/";
    my @list = ('me2s0acc1','me2s0acc2','me2s0acc3','me2s0acc4',
		'me2s0acc5','me2s0acc6','me2s0acc7','me2s0acc8','me2s0acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	#print STDERR $file."\t".$num."\t".$n."\n";
	$num++;
    }
    return @metables;
}
sub maxentscore{
    my $seq = shift;
    my $table_ref = shift;
    my @metables = @$table_ref;
    my @sc;
    $sc[0] = $metables[0]{&hashseq(substr($seq,0,7))};
    $sc[1] = $metables[1]{&hashseq(substr($seq,7,7))};
    $sc[2] = $metables[2]{&hashseq(substr($seq,14,7))};
    $sc[3] = $metables[3]{&hashseq(substr($seq,4,7))};
    $sc[4] = $metables[4]{&hashseq(substr($seq,11,7))};
    $sc[5] = $metables[5]{&hashseq(substr($seq,4,3))};
    $sc[6] = $metables[6]{&hashseq(substr($seq,7,4))};
    $sc[7] = $metables[7]{&hashseq(substr($seq,11,3))};
    $sc[8] = $metables[8]{&hashseq(substr($seq,14,4))};
    my $finalscore = $sc[0] * $sc[1] * $sc[2] * $sc[3] * $sc[4] / ($sc[5] * $sc[6] * $sc[7] * $sc[8]);
    return $finalscore;
}    
    


sub getrest{
  my $seq = shift;
  my $seq_noconsensus = substr($seq,0,18).substr($seq,20,3);
  return $seq_noconsensus;
}

sub scoreconsensus{
  my $seq = shift;
  my @seqa = split(//,uc($seq));
  my %bgd; 
  $bgd{'A'} = 0.27; 
  $bgd{'C'} = 0.23; 
  $bgd{'G'} = 0.23; 
  $bgd{'T'} = 0.27;  
  my %cons1;
  $cons1{'A'} = 0.9903;
  $cons1{'C'} = 0.0032;
  $cons1{'G'} = 0.0034;
  $cons1{'T'} = 0.0030;
  my %cons2;
  $cons2{'A'} = 0.0027; 
  $cons2{'C'} = 0.0037; 
  $cons2{'G'} = 0.9905; 
  $cons2{'T'} = 0.0030;
  my $addscore = $cons1{$seqa[18]} * $cons2{$seqa[19]}/ ($bgd{$seqa[18]} * $bgd{$seqa[19]}); 
  return $addscore;
}
sub log2{
      my ($val) = @_;
    return log($val)/log(2);
}
