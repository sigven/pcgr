use strict;


my $inputfile = $ARGV[0];
my $usemaxent = 1;

my %me2x5 = &makescorematrix('me2x5');
my %seq = &makesequencematrix('splicemodels/splice5sequences');

my %bgd;
$bgd{'A'} = 0.27;
$bgd{'C'} = 0.23;
$bgd{'G'} = 0.23;
$bgd{'T'} = 0.27; 



open (FILE,"<$inputfile") || die "can't open!\n";

while(<FILE>) {
    chomp;
    if (/^\s*$/) { #discard blank lines;
	next;
    } 
    elsif (/^>/) { #discard comment lines;
	next;
    }
    else {
        $_ =~ s/\cM//g; #gets rid of carriage return
	my $str = $_;
	print $str."\t";
	$str = uc($str);
	if ($usemaxent) { 
	 print sprintf("%.2f",&log2(&scoreconsensus($str)*$me2x5{$seq{&getrest($str)}}))."\n";
	}
    }
}

  
sub makesequencematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$_} = $n;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}
sub makescorematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$n} = $_;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}

sub getrest{
  my $seq = shift;
  my @seqa = split(//,uc($seq));
  return $seqa[0].$seqa[1].$seqa[2].$seqa[5].$seqa[6].$seqa[7].$seqa[8];
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
  $cons1{'A'} = 0.004;
  $cons1{'C'} = 0.0032;
  $cons1{'G'} = 0.9896;
  $cons1{'T'} = 0.0032;
  my %cons2;
  $cons2{'A'} = 0.0034; 
  $cons2{'C'} = 0.0039; 
  $cons2{'G'} = 0.0042; 
  $cons2{'T'} = 0.9884;
  my $addscore = $cons1{$seqa[3]}*$cons2{$seqa[4]}/($bgd{$seqa[3]}*$bgd{$seqa[4]}); 
  return $addscore;
}

sub log2{
      my ($val) = @_;
    return log($val)/log(2);
}
