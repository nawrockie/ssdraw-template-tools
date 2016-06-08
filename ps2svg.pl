#!/usr/bin/env perl

use strict;
use Getopt::Long;

&GetOptions();

my $usage;
$usage  = "ps2svg.pl [OPTIONS] <esl-ssdraw postscript template file>\n";
if(scalar(@ARGV) != 1) { die $usage; }

my $psfile = $ARGV[0];

my $psIFH;
open($psIFH, $psfile) || die "ERROR unable to open $psfile";

# globals
our $fontG = "Courier";
our $fontsizeG = "8";

# nucleotide data
my @nt_nA = (); # actual nucleotide
my @nt_xA = (); # x coord
my @nt_yA = (); # y coord

# bpconnect
my @bp_x1A = (); # first x coord
my @bp_x2A = (); # second x coord
my @bp_y1A = (); # first y coord
my @bp_y2A = (); # second y coord

# Input
my $line;
my $scale = "";
while($line = <$psIFH>) { 
  chomp $line;
  if($line =~ /^(\d+\.*\d*)\s+(\d+\.*\d*)\s+scale$/) { 
    $scale = $1; 
    if($scale ne $2) { die "ERROR illegal scale line: $line"; }
  }
  if($line =~ m/\% begin text nucleotides/) { 
    ps_input_nucleotides($psIFH, \@nt_nA, \@nt_xA, \@nt_yA);
  }
  if($line =~ m/\% begin lines bpconnects/) { 
    ps_input_bpconnects($psIFH, \@bp_x1A, \@bp_y1A, \@bp_x2A, \@bp_y2A);
  }
}
close($psIFH);
if($scale eq "") { die "ERROR did not read postscript scale line\n"; }
# Output 

svg_output_header($scale);
svg_output_nucleotides($scale, \@nt_nA, \@nt_xA, \@nt_yA);
svg_output_bpconnects($scale, \@bp_x1A, \@bp_y1A, \@bp_x2A, \@bp_y2A);
svg_output_tail();

# Subroutines:

sub ps_input_nucleotides {
  my ($psIFH, $nt_nAR, $nt_xAR, $nt_yAR) = @_;
  my $i = 0;
  @{$nt_nAR} = ();
  @{$nt_xAR} = ();
  @{$nt_yAR} = ();
  while($line = <$psIFH>) { 
    # (G) 168.00 392.00 moveto show
    if($line =~ /^\((\S)\)\s+(\-?\d+\.*\d*)\s+(\-?\d+\.*\d*)\s+moveto\s+show/) { 
      ($nt_nAR->[$i], $nt_xAR->[$i], $nt_yAR->[$i]) = ($1, $2, $3);
      $i++;
    }
    if($line =~ m/\% end text nucleotides/) { last; }
  }
}

sub ps_input_bpconnects {
  my ($psIFH, $bp_x1AR, $bp_y1AR, $bp_x2AR, $bp_y2AR) = @_;
  my $i = 0;
  @{$bp_x1AR} = ();
  @{$bp_x2AR} = ();
  @{$bp_y1AR} = ();
  @{$bp_y2AR} = ();
  while($line = <$psIFH>) { 
    # 179.00 394.00 175.00 394.00 newpath moveto lineto stroke
    if($line =~ /^(\-?\d+\.*\d*)\s+(\-?\d+\.*\d*)\s+(\-?\d+\.*\d*)\s+(\-?\d+\.*\d*)\s+newpath moveto lineto stroke/) { 
      ($bp_x1AR->[$i], $bp_y1AR->[$i], $bp_x2AR->[$i], $bp_y2AR->[$i]) = ($1, $2, $3, $4);
      $i++;
    }
    if($line =~ m/\% end lines bpconnects/) { last; }
  }
}

sub svg_output_nucleotides {
  my ($scale, $nt_nAR, $nt_xAR, $nt_yAR) = @_;
  my $i = 0;
  my $nnt = scalar(@{$nt_nAR});
  my $height = 792;
  if($nnt != scalar(@{$nt_xAR})) { die "ERROR nt_xAR has different size than nt_ntAR"; }
  if($nnt != scalar(@{$nt_yAR})) { die "ERROR nt_yAR has different size than nt_ntAR"; }

  for ($i = 0; $i < $nnt; $i++) {
    svg_output_single_nucleotide($i, $scale * $nt_xAR->[$i], $height - ($scale * $nt_yAR->[$i]), $nt_nAR->[$i], $fontsizeG * $scale);
  }
}

sub svg_output_single_nucleotide {
  my ($i, $x, $y, $nt, $fontsize) = @_;

  print ("<text");
  printf("\n\tid=\"%s\"", "nt" . ($i+1));
  printf("\n\tx=\"%s\"", $x);
  printf("\n\ty=\"%s\"", $y);
  printf("\n\tstyle=\"font-size:%spx;font-family:$fontG\"", $fontsize);
  printf(">%s<", $nt);
  print ("/text>\n");
}

sub svg_output_bpconnects {
  my ($scale, $bp_x1AR, $bp_y1AR, $bp_x2AR, $bp_y2AR) = @_;
  my $i = 0;
  my $nbp = scalar(@{$bp_x1AR});
  my $height = 792;
  if($nbp != scalar(@{$bp_x2AR})) { die "ERROR bp_x2AR has different size than bp_x1AR"; }
  if($nbp != scalar(@{$bp_y1AR})) { die "ERROR bp_y1AR has different size than bp_x1AR"; }
  if($nbp != scalar(@{$bp_y2AR})) { die "ERROR bp_y2AR has different size than bp_x1AR"; }

  for ($i = 0; $i < $nbp; $i++) {
    svg_output_single_bpconnect($i, $scale * $bp_x1AR->[$i], $height - ($scale * $bp_y1AR->[$i]), $scale * $bp_x2AR->[$i], $height - ($scale * $bp_y2AR->[$i]));
  }
}

sub svg_output_single_bpconnect {
  my ($i, $x1, $y1, $x2, $y2) = @_;

  print ("<path");
  printf("\n\tstyle=\"fill:none;stroke:\#000000;stroke-width:1px;stroke-opacity:1\"");
  printf("\n\td=\"M %.2f %.2f L %.2f %.2f\"", $x1, $y1, $x2, $y2);
  printf("\n\tid=\"%s\"", "bp" . ($i+1));
  printf("/>\n");
}

sub svg_output_header { 
  print("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
  print("<!-- Generator: Adobe Illustrator 14.0.0, SVG Export Plug-In . SVG Version: 6.00 Build 43363)  -->\n");
  print("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
  print("<svg version=\"1.1\" id=\"Layer_1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" x=\"0px\" y=\"0px\"\n");
  print("width=\"612px\" height=\"792px\" viewBox=\"0 0 612 792\" enable-background=\"new 0 0 612 792\" xml:space=\"preserve\" scale=\"$scale\" \>\n");
}

sub svg_output_tail {
  print("</svg>\n");
}
