=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
    
=cut

=head1 NAME

 NearestExonJB

=head1 SYNOPSIS

 mv NearestExonJB.pm ~/.vep/Plugins
 ./vep -i variations.vcf --cache --plugin NearestExonJB

=head1 DESCRIPTION

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that
 finds the nearest exon junction boundary to a coding sequence variant. More than 
 one boundary may be reported if the boundaries are equidistant.

 The plugin will report the Ensembl identifier of the exon, the distance to the
 exon boundary, the boundary type (start or end of exon) and the total
 length in nucleotides of the exon.

 Various parameters can be altered by passing them to the plugin command:

 - max_range : maximum search range in bp (default: 10000)

 Parameters are passed e.g.:

 --plugin NearestExonJB,max_range=50000

=cut

package NearestExonJB;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

my $char_sep = "|";

my %CONFIG = (
  max_range => 10000,
);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my $params = $self->params;

  # get output format
  $char_sep = "+" if ($self->{config}->{output_format} eq 'vcf');

  foreach my $param(@$params) {
    my ($key, $val) = split('=', $param);
    die("ERROR: Failed to parse parameter $param\n") unless defined($key) && defined($val);
    $CONFIG{$key} = $val;
  }

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $header = 'Nearest Exon Junction Boundary (coding sequence variants only). Format:';
  $header .= join($char_sep, qw(ExonID distance start/end length) );

  return {
    NearestExonJB => $header,
  }
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->base_variation_feature;
  my $trv = $tva->base_transcript_variation;

  my $loc_string = sprintf("%s:%s-%i-%i", $trv->transcript_stable_id, $vf->{chr} || $vf->seq_region_name, $vf->{start}, $vf->{end});

  if(!exists($self->{_cache}) || !exists($self->{_cache}->{$loc_string})) {
    my $exons = $trv->_overlapped_exons;
    my %dists;
    my $min = $CONFIG{max_range};
    foreach my $exon (@$exons) {
      my $startD = abs ($vf->start - $exon->seq_region_start);
      my $endD = abs ($vf->end - $exon->seq_region_end);
      if ($startD < $endD){
        $dists{$exon->stable_id}{$startD} = 'start';
        $dists{$exon->stable_id}{len} = $exon->length;
        $min = $startD if $min > $startD;
      } elsif ($startD > $endD){
        $dists{$exon->stable_id}{$endD} = 'end';
        $dists{$exon->stable_id}{len} = $exon->length;
        $min = $endD if $min > $endD;
      } else {
        $dists{$exon->stable_id}{$startD} = "start_end";
        $dists{$exon->stable_id}{len} = $exon->length;
        $min = $startD if $min > $startD;
      }
    }

    my @finalRes;
    foreach my $exon (keys %dists){
      if (exists $dists{$exon}{$min}) {
        push(@finalRes, $exon.$char_sep.$min.$char_sep.$dists{$exon}{$min}.$char_sep.$dists{$exon}{len})
      }
    }

    $self->{_cache}->{$loc_string} = scalar @finalRes ? join(",", @finalRes) : undef;
  }
  return $self->{_cache}->{$loc_string} ? { NearestExonJB => $self->{_cache}->{$loc_string} } : {};
}

1;

