#!/usr/bin/env perl
# 
# EPN, Fri Apr 12 13:44:52 2013
# svg2ps.pl
# 
# Convert a SVG file created by ps2svg.pl and probably modified in 
# Inkscape (by only duplicating existing nucleotides and moving them
# or by deleting nucleotides) back into postscript for esl-ssdraw.
# 
# Based partly on ssd_addstruct.pl [EPN, Sun Jun  7 14:35:56 2009]
#

use strict;
use Getopt::Long;
use Math::Trig;

my $do_leave_seq = 0;
my $do_overwrite_seq = 1;
my $do_rf_seq = 0;
my $input_is_svg = 1;
my $input_is_crw = 0;
my $textmod = "";

&GetOptions( "L"   => \$do_leave_seq, 
             "R"   => \$do_rf_seq,
             "C"   => \$input_is_crw,
             "T=s" => \$textmod);
if($do_leave_seq && $do_rf_seq) { die "ERROR only one of -L and -R can be used"; }
if($do_leave_seq) { 
  $do_overwrite_seq = 0; 
}
if($input_is_crw) { $input_is_svg = 0; }

our $fontsize = 8;

my $usage;
$usage  = "svg2ps.pl [OPTIONS] <svg file> <original esl-ssdraw postscript template file> <Pfam format SS_cons annotated alignment file>\n\t";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t-L:     leave SS diagram sequence alone; don't change it\n";
$usage .= "\t\t-R:     overwrite SS diagram sequence with RF sequence from alignment\n";
$usage .= "\t\t-C:     <svg file> is actually a CRW postscript file, NOT an SVG file\n";
$usage .= "\t\t-T <n>: label every <n>'th nucleotide [default: autodetermine]\n";
$usage .= "\n\n";

if(scalar(@ARGV) != 3) { die $usage; }

my ($svg_or_crwfile, $orig_psfile, $stkfile) = @ARGV;

my $svg_or_crwIFH;
open($svg_or_crwIFH, $svg_or_crwfile) || die "ERROR unable to open $svg_or_crwfile";

my $psIFH;
open($psIFH, $orig_psfile) || die "ERROR unable to open $orig_psfile";

my $stkIFH;
open($stkIFH, $stkfile) || die "ERROR unable to open $stkfile";

# get original ps file data
# we will regurgitate the entire file, but replace the following sections:
# 'text positiontext'
# 'lines positionticks'
# 'lines bpconnects'
# 'text nucleotides'
my @orig_pslinesA = ();
my $orig_ps_line_postext = 0;   # line in @orig_pslinesA after which 'text positiontext' section appeared
my $orig_ps_line_posticks = 0;  # line in @orig_pslinesA after which 'lines positionticks' section appeared
my $orig_ps_line_bps = 0;       # line in @orig_pslinesA after which 'lines bpconnects' section appeared
my $orig_ps_line_nts = 0;       # line in @orig_pslinesA after which 'text nucleotides' section appeared
ps_input_original($psIFH, \@orig_pslinesA, \$orig_ps_line_postext, \$orig_ps_line_posticks, \$orig_ps_line_bps, \$orig_ps_line_nts);
close($psIFH);

# get basepair info from stockholm alignment
my @bp_ctA = ();
my $clen = 0;
my ($seq, $rfseq);
stk_input($stkIFH, \@bp_ctA, \$seq, \$rfseq, \$clen);
if($do_rf_seq && $rfseq eq "") { die "ERROR with -R, alignment must have GC RF annotation\n"; }

my @nt_nA = (); # actual nucleotide
my @nt_xA = (); # x coord
my @nt_yA = (); # y coord

# input either SVG file (default) or CRW postscript file (if -C enabled)
if($input_is_svg){ 
  input_svg_file($svg_or_crwIFH, $clen, \@nt_nA, \@nt_xA, \@nt_yA);
}
else { # read CRW postscript file
  input_crw_ps_file($svg_or_crwIFH, $clen, \@nt_nA, \@nt_xA, \@nt_yA);
}

# overwrite @nt_nA if necessary
if($do_overwrite_seq) { 
  @nt_nA = split("", $seq);
}
if($do_rf_seq) { 
  @nt_nA = split("", $rfseq);
}

# pre-generate output lines for base pairs and tick marks
my @bp_xbA = (); # array of basepairs begin x coordinate
my @bp_xeA = (); # array of basepairs end   x coordinate
my @bp_ybA = (); # array of basepairs begin y coordinate
my @bp_yeA = (); # array of basepairs end   y coordinate
my @bp_toprintA = ();   # lines to output in bpconnects section
my @postick_toprintA = (); # lines to output in positionticks section
my @postext_toprintA = (); # '100' text output in positiontext section

my ($nt_i, $nt_j);
if($textmod eq "") { 
  $textmod = 10;
  if($clen >= 60)  { $textmod = 20; }
  if($clen >= 150) { $textmod = 50; }
  if($clen >= 500) { $textmod = 100; }
}

# first, determine bp info and output lines
my $nbp = 0;
my ($i, $j, $line);
for($i = 1; $i <= $clen; $i++) { 
  $nt_i = ($i-1);
  $j = $bp_ctA[$i];
  if($j > $i) { 
    $nt_j = $j-1;
    print_bp_line(\@bp_toprintA, \@bp_xbA, \@bp_xeA, \@bp_ybA, \@bp_yeA, $nt_xA[$nt_j], $nt_yA[$nt_j], $nt_xA[$nt_i], $nt_yA[$nt_i]);
    $nbp++;
  }
}

# second pass through, determine tick and position text output lines
for($i = 1; $i <= $clen; $i++) { 
  if($i % 10 == 0) { 
    print_tick_number(\@postick_toprintA, \@postext_toprintA, $textmod, $i-1, \@nt_xA, \@nt_yA, \@bp_xbA, \@bp_xeA, \@bp_ybA, \@bp_yeA);
  }
}

# Output 
my $orig_nlines = scalar(@orig_pslinesA);
for($i = 0; $i < $orig_nlines; $i++) { 
  print $orig_pslinesA[$i];
  if($i == $orig_ps_line_postext)  { ps_output_postext    (\@postext_toprintA); }
  if($i == $orig_ps_line_posticks) { ps_output_posticks   (\@postick_toprintA); }
  if($i == $orig_ps_line_bps)      { ps_output_bpconnects (\@bp_toprintA); }
  if($i == $orig_ps_line_nts)      { ps_output_nucleotides(\@nt_nA, \@nt_xA, \@nt_yA); }
}

# Subroutines:

sub ps_input_original { 
  my ($psIFH, $linesAR, $line_postextR, $line_posticksR, $line_bpsR, $line_ntsR)= @_;

  my $ignore_flag = 0;
  my $line_ctr = 0;
  my $nfound = 0; # this should be 8 by end of function
  while($line = <$psIFH>) { 
    if($line =~ m/^% begin text positiontext/) { 
      $ignore_flag = 1; 
      $$line_postextR = $line_ctr;
      $nfound++;
    }
    if($line =~ m/^% begin lines positionticks/) { 
      $ignore_flag = 1; 
      $$line_posticksR = $line_ctr;
      $nfound++;
    }
    if($line =~ m/^% begin lines bpconnects/) { 
      $ignore_flag = 1; 
      $$line_bpsR = $line_ctr;
      $nfound++;
    }
    if($line =~ m/^% begin text nucleotides/) { 
      $ignore_flag = 1; 
      $$line_ntsR = $line_ctr;
      $nfound++;
    }

    # store line, if ignore_flag is down
    if(! $ignore_flag) { 
      $linesAR->[$line_ctr] = $line;
      $line_ctr++;
    }

    if($line =~ m/^% end text positiontext/) { 
      $ignore_flag = 0; 
      $nfound++;
    }
    if($line =~ m/^% end lines positionticks/) { 
      $ignore_flag = 0; 
      $nfound++;
    }
    if($line =~ m/^% end lines bpconnects/) { 
      $ignore_flag = 0; 
      $nfound++;
    }
    if($line =~ m/^% end text nucleotides/) { 
      $ignore_flag = 0; 
      $nfound++;
    }
  }
  if($nfound != 8) { die "ERROR didn't find all 8 special lines in ps_parse_original(), something wrong with PS file?"; }
}

sub svg_input_scale { 
  my ($svgIFH, $line) = @_;

#	<svg
#        OTHER CRAP
#        scale="1.7"
#        OTHER CRAP
#       />

  my $keep_going = 1;
  my $scale = "";
  while($keep_going) { 
    if($line =~ m/scale=\"(\d+\.*\d*)\"/) { $scale = $1; $keep_going = 0; }
    if($line =~ m/\/>/)                   { $keep_going = 0; }
    if($keep_going) { 
      if(! ($line = <$svgIFH>)) { die "ERROR ran out of file when parsing SVG header"; }
    }
  }
  if($scale eq "") { die "ERROR didn't read scale in SVG header"; }
  return $scale;
}

sub svg_input_nucleotide { 
  my ($svgIFH, $ntR, $xR, $yR, $idR) = @_;
#
# If nt is unmodified from ps2svg.pl conversion:
#	id="nt4"
#	x="285.6"
#	y="166.4"
#	style="font-size:13.6px;font-family:Helvetica">G</text>
#
# If nt was created by duplicating an existing one in Inkscape
#    style="font-size:13.60000038px;font-family:Helvetica"
#       y="248"
#       x="156.3"
#      id="nt17.1">U</text>
#
# This function is a bit more permissive, it allows any 
# order for the 4 lines and for the '>U</text>' to be added at the end of any
  my $keep_going;
  my $seen_id    = 0;
  my $seen_x     = 0;
  my $seen_y     = 0;
  my $seen_style = 0;
  my $seen_nt    = 0;

  while((! $seen_nt) && ($line = <$svgIFH>)) { 
    if   ($line =~ /^\s*id=\"nt(\d+\.*\d*)\"/) { $$idR = $1; $seen_id = 1; }
    elsif($line =~ /^\s*x=\"(\d+\.*\d*)\"/)    { $$xR = $1;  $seen_x  = 1; }
    elsif($line =~ /^\s*y=\"(\d+\.*\d*)\"/)    { $$yR = $1;  $seen_y  = 1; }
    elsif($line =~ /^\s*style/)                {             $seen_style = 1; } # we don't parse style line
    elsif($line =~ /^\s*transform/)            { ; } # we don't parse transform line
    else { die "ERROR did not parse SVG nucleotide line $line"; }

    # nucleotide is appended to end of one of the four lines
    if($line =~ />(\S)\<\/text>/) { $$ntR = $1; $seen_nt = 1; }
  }

  if(! $seen_id) { die "ERROR failed to parse id in SVG nucleotide"; }
  if(! $seen_x)  { die "ERROR failed to parse id in SVG x"; }
  if(! $seen_y)  { die "ERROR failed to parse id in SVG y"; }
  if(! $seen_nt) { die "ERROR failed to parse id in SVG nt"; }
  #printf("$$idR $$xR $$yR $$ntR\n");

  return;
}

sub stk_input { 
  my ($stkIFH, $bp_ctAR, $seqR, $rfseqR, $clenR) = @_;

  my $seen_struct = 0;
  my $seen_seq    = 0;
  my $seen_rf = 0;
  my $firstseq = "";
  my $seqname = "";
  my $ss = "";
  my $rf = "";
  while($line = <$stkIFH>) { 
    if(! $seen_seq) { 
      if($line =~ m/\w/ && $line =~ m/^[^\#]/) { 
        if($line =~ /(\S+)\s+(\S+)/) { 
          $seqname = $1;
          $firstseq = $2;
          $seen_seq = 1;
        }
        else { 
          die "ERROR unexpected line: $line in stockholm alignment"; 
        }
      }
    }
    else { # $seen_seq is true
      if($line =~ /^($seqname)\s+(\S+)/) { 
        $firstseq .= $2;
      }
    }
    if($line =~ s/^\#=GC SS\_cons\s+//) {
      chomp $line; 
      $ss .= $line;
      $seen_struct = 1;
    }
    if($line =~ s/^\#=GC RF\s+//) {
      chomp $line; 
      $rf .= $line;
      $seen_rf = 1;
    }
  }
  close($stkIFH);
  if(! $seen_struct) { 
    die "ERROR did not read SS_cons information"; 
  }

  struct_seq2bp_arr($ss, $bp_ctAR, $clenR);

  $$seqR = $firstseq;
  $$rfseqR = $rf;
  return;
}

sub ps_output_nucleotides {
  my ($nt_nAR, $nt_xAR, $nt_yAR) = @_;
  my $i = 0;
  my $nnt = scalar(@{$nt_nAR});
  if($nnt != scalar(@{$nt_xAR})) { die "ERROR nt_xAR has different size than nt_ntAR"; }
  if($nnt != scalar(@{$nt_yAR})) { die "ERROR nt_yAR has different size than nt_ntAR"; }

  print "% begin text nucleotides\n";

  for ($i = 0; $i < $nnt; $i++) {
    printf("(%s) %.2f %.2f moveto show\n",
           $nt_nAR->[$i], $nt_xAR->[$i], $nt_yAR->[$i]);
    if((($i+2) % 10 == 0) && ($i < ($nnt-1))) { 
      printf("%% nucleotide %d on next line\n", ($i+2)); 
    }
  }
  print "% end text nucleotides\n";
}

sub ps_output_bpconnects { 
  my ($toprintAR) = $_[0];

  print "% begin lines bpconnects\n";
  foreach my $line (@{$toprintAR}) { 
    print $line; 
  }
  print "% end lines bpconnects\n";
}

sub ps_output_posticks { 
  my ($toprintAR) = $_[0];

  print "% begin lines positionticks\n";
  foreach my $line (@{$toprintAR}) { 
    print $line; 
  }
  print "% end lines positionticks\n";
}

sub ps_output_postext { 
  my ($toprintAR) = $_[0];

  print "% begin text positiontext\n";
  foreach my $line (@{$toprintAR}) { 
    print $line; 
  }
  print "% end text positiontext\n";
}

#################################################################
# subroutine : struct_seq2bp_arr
# sub class  : sequence
#
# EPN 05.02.06
# 
# purpose : Given a structure sequence (where "." are treated
#           as gaps), fill a bp array. Print an error and 
#           die if the structure sequence is not well-nested.
#           sequence that has a gap for either i or j (or both).
#
# args : (1) $struct_seq
#            the structure sequence ("."'s read as gaps)
#        (2) $bp_AR
#            array to the bp array to fill
#        (3) $clenR
#            ref to consensus length, to fill
#################################################################
sub struct_seq2bp_arr
{
  my ($struct_seq, $bp_AR, $clenR) = @_;
  
  #left to '<'
  $struct_seq=~ s/[\(\[\{]/</g;
  #right to '>'
  $struct_seq=~ s/[\)\]\}]/>/g;
  #single stranded to :
  $struct_seq=~ s/[\,\-\_]/\:/g;
  #everything else to a '.'
  $struct_seq=~ s/[^\.><\:]/\./g;
  
  my @struct_A = split("", ("!" . $struct_seq));
  
  my $clen = scalar(@struct_A) - 1;
  
  my @stack = ();
  my ($i, $left, $right); 
  
  $bp_AR->[0] = -1;
  for($i = 1; $i < $clen; $i++) { 
    $bp_AR->[$i] = -1;
  }
  for($i = 1; $i < $clen; $i++) {
    if($struct_A[$i] ne "\.") { 
      if($struct_A[$i] eq "<") { 
        #left side
        push(@stack, $i);
      }
      elsif($struct_A[$i] eq ">") { 
        #right side
        if(scalar(@stack) < 1) { 
          #trying to pop empty stack
          die "ERROR : trying to pop empty stack for i: $i";
        }
        $left = pop @stack;
        $right = $i;
        $bp_AR->[$left] = $right;
        $bp_AR->[$right] = $left;
      }
    }
    else { #single stranded  
      $bp_AR->[$i] = 0;
    }
  }
  $$clenR = $clen;
  return;
}

############################################################
sub print_bp_line
{
  my ($bp_toprintAR, $bp_xbAR, $bp_xeAR, $bp_ybAR, $bp_yeAR, $xi, $yi, $xj, $yj) = @_;
  
  my $fontsize = 8;
  my $bpscalar = 1.0;
  my $hw = 0.75;

  my $xdiff = abs($xi - $xj);
  my $ydiff = abs($yi - $yj);
  my $xeq = ($xdiff < 0.0001) ? 1 : 0; 

  # degree of angle separating i and j residues
  my $deg = ($xeq) ? 90. : rad2deg(atan(($ydiff/$xdiff)));

  # find middle of x and y residues
  my $mid_xi = $xi + ($fontsize * $hw / 2.);
  my $mid_xj = $xj + ($fontsize * $hw / 2.);
  my $mid_yi = $yi + ($fontsize * $hw / 2.);
  my $mid_yj = $yj + ($fontsize * $hw / 2.);

  # find distance to skip, so we don't overlap the residue
  my $xskip = $bpscalar * $fontsize * cos(deg2rad($deg));
  my $yskip = $bpscalar * $fontsize * sin(deg2rad($deg));

  # finally, get begin and end coords for x and y, by subtracting/adding {x,y}skip
  my ($xb, $xe, $yb, $ye);
  if($xi > $xj) { 
    $xb = $mid_xi - $xskip;
    $xe = $mid_xj + $xskip;
  }
  else { 
    $xb = $mid_xi + $xskip;
    $xe = $mid_xj - $xskip;
  }
  if($yi > $yj) { 
    $yb = $mid_yi - $yskip;
    $ye = $mid_yj + $yskip;
  }
  else { 
    $yb = $mid_yi + $yskip;
    $ye = $mid_yj - $yskip;
  }
  
  push(@{$bp_toprintAR}, sprintf("%.2f %.2f %.2f %.2f newpath moveto lineto stroke\n", $xb, $yb, $xe, $ye));
  #push(@{$bp_toprintAR}, sprintf("%.2f %.2f %.2f %.2f newpath moveto lineto stroke\n", $mid_xi, $mid_yi, $mid_xj, $mid_yj));
  push(@{$bp_xbAR}, $xb);
  push(@{$bp_xeAR}, $xe);
  push(@{$bp_ybAR}, $yb);
  push(@{$bp_yeAR}, $ye);

  return;
}

############################################################
sub print_tick_number
{
  my ($postick_toprintAR, $postext_toprintAR, $textmod, $num, $xAR, $yAR, $bpxbAR, $bpxeAR, $bpybAR, $bpyeAR) = @_;
  my $clen = scalar(@{$xAR});
  my ($x, $y, $xprv, $yprv, $xnxt, $ynxt);
  my $i = $num;
  my $num2print = $num+1; # off by-one in array indexing
  
  my $space = 3.;
  my $mult  = 3.5;
  my $xprv = $xAR->[($i-1)];
  my $yprv = $yAR->[($i-1)];
  my $x    = $xAR->[$i];
  my $y    = $yAR->[$i];
  if($i == ($clen-1)) { $xnxt = $x; $ynxt = $y; }
  else { 
    $xnxt = $xAR->[($i+1)];
    $ynxt = $yAR->[($i+1)];
  }
  #printf("(%s) %.2f %.2f lwstring\n", $num, $x+5., $y+5.);

  # get angle by looking at nearest neighboring residue on each side
  if(($xnxt eq "end") && ($ynxt eq "end")) { $xnxt = $x; $ynxt = $y; }
  my $deg = angle_between_two_points($xnxt, $ynxt, $xprv, $yprv);
  my $sin = sin(deg2rad($deg));
  my $cos = cos(deg2rad($deg));
  my $d1 = ($space) * $sin;
  my $d2 = ($space) * $cos;

  #printf STDERR ("TMP %4d %.2f %.2f deg: %.2f d1: %.2f d2: %.2f\n",  $num, $x+($fontsize/2.), $y+($fontsize/2.), $deg, $d1, $d2);

  my ($xa1, $xa2, $xa, $xb1, $xb2, $xb);
  my ($ya1, $ya2, $ya, $yb1, $yb2, $yb);
  my $xinc = 0;
  my $yinc = 0;
  if($xnxt > $xprv) { $xinc = 1; }
  if($ynxt > $yprv) { $yinc = 1; }
  if($xinc) { 
    $xa1 = $x+($space) - $d1;
    $xa2 = $x+($space) - (($mult)*$d1);
    $xa  = $x+($space) - ((($mult+1)/2.) * $d1);
    $xb1 = $x+($space) + $d1;
    $xb2 = $x+($space) + (($mult)*$d1);
    $xb  = $x+($space) + ((($mult+1)/2.) * $d1);
  }
  else { 
    $xa1 = $x+($space) + $d1;
    $xa2 = $x+($space) + (($mult)*$d1);
    $xa  = $x+($space) + ((($mult+1)/2.) * $d1);
    $xb1 = $x+($space) - $d1;
    $xb2 = $x+($space) - (($mult)*$d1);
    $xb  = $x+($space) - ((($mult+1)/2.) * $d1);
  }	
  if($yinc) { 
    $ya1 = $y+($space) + $d2;
    $ya2 = $y+($space) + (($mult)*$d2);
    $ya  = $y+($space) + ((($mult+1)/2.) * $d2);
    $yb1 = $y+($space) - $d2;
    $yb2 = $y+($space) - (($mult)*$d2);
    $yb  = $y+($space) - ((($mult+1)/2.) * $d2);
  }
  else { 
    $ya1 = $y+($space) - $d2;
    $ya2 = $y+($space) - (($mult)*$d2);
    $ya  = $y+($space) - ((($mult+1)/2.) * $d2);
    $yb1 = $y+($space) + $d2;
    $yb2 = $y+($space) + (($mult)*$d2);
    $yb  = $y+($space) + ((($mult+1)/2.) * $d2);
  }
  my $resdista = dist2closest_other_residue($xAR, $yAR, $i+1, $xa, $ya);
  my $resdistb = dist2closest_other_residue($xAR, $yAR, $i+1, $xb, $yb);
  my $bpdista  = dist2closest_basepair($bpxbAR, $bpxeAR, $bpybAR, $bpyeAR, $xa, $ya);
  my $bpdistb  = dist2closest_basepair($bpxbAR, $bpxeAR, $bpybAR, $bpyeAR, $xb, $yb);

  my $dista = $resdista;
  if($bpdista < $dista) { $dista = $bpdista; } 
  my $distb = $resdistb;
  if($bpdistb < $distb) { $distb = $bpdistb; } 

  #printf STDERR ("TMP %4d $resdista $bpdista    $resdistb $bpdistb      $dista  $distb\n\n", $i);

  my ($xnum, $ynum);
  if($dista > $distb) { # a is furthest from it's closest other residue
    push(@{$postick_toprintAR}, sprintf("%.2f %.2f %.2f %.2f newpath moveto lineto stroke\n", $xa1, $ya1, $xa2, $ya2));
    $xnum = $xa2;
    $ynum = $ya2;
  }
  else { # b is furthest from it's closest other residue
    push(@{$postick_toprintAR}, sprintf("%.2f %.2f %.2f %.2f newpath moveto lineto stroke\n", $xb1, $yb1, $xb2, $yb2));
    $xnum = $xb2;
    $ynum = $yb2;
  }

  my $tmp;
  my $ndig;
  my ($xeq, $yeq);
  if(($i+1) % $textmod == 0) { # draw a number of this position next to the tick mark
    # There's 8 different scenarios we consider:
    # 
    # The 8 positions:                       Example of where number goes:
    #                     8  1  2                         50  50 50
    #                      \ | /                            \ | /
    #                       \|/                              \|/
    #                   7---- ----3                     50---- ----50
    #                       /|\                              /|\
    #                      / | \                            / | \
    #                     6  5  4                         50  50 50
    # 
    # Values per position:
    #             1    2    3    4    5    6    7    8
    # xeq:        yes  no   no   no   yes  no   yes  no          
    # yeq:        no   no   yes  no   no   no   yes  no
    # xinc:       yes  yes  yes  yes  yes  no   no   no
    # yinc:       yes  yes  yes  no   no   no   yes  yes
    #
    if($dista > $distb) { # a is furthest from it's closest other residue
      $xnum = $xa2;
      $ynum = $ya2;
      $xinc = ($xa2 >= $xa1) ? 1 : 0;
      $yinc = ($ya2 >= $ya1) ? 1 : 0;
      $xeq  = (abs($xa2 - $xa1) < 0.0001) ? 1 : 0;
      $yeq  = (abs($ya2 - $ya1) < 0.0001) ? 1 : 0;
    }
    else { 
      $xnum = $xb2;
      $ynum = $yb2;
      $xinc = ($xb2 >= $xb1) ? 1 : 0;
      $yinc = ($yb2 >= $yb1) ? 1 : 0;
      $xeq  = (abs($xb2 - $xb1) < 0.0001) ? 1 : 0;
      $yeq  = (abs($yb2 - $yb1) < 0.0001) ? 1 : 0;
    }
    $ndig = 1;
    $tmp  = $num2print;
    while($tmp > 10) { $tmp/=10; $ndig++; }
    if($xeq)       { $xnum -= 0.55 * $fontsize; }
    elsif($xinc)   { $xnum += 0.05  * $fontsize; }
    else           { $xnum -= 0.55 * $ndig * $fontsize; } # alternatively: add '* $cos'

    if($yeq)       { $ynum -= 0.375 * $fontsize; }
    elsif($yinc)   { $ynum += 0.05 * $fontsize; }
    else           { $ynum -= 0.75 * $fontsize; }
    push(@{$postext_toprintAR}, sprintf("(%s) %.2f %.2f moveto show\n", $num2print, $xnum, $ynum));
  }
}
############################################################
sub dist2closest_basepair
{
  my ($bpxbAR, $bpxeAR, $bpybAR, $bpyeAR, $x, $y) = @_;
  my $mindist = 100000;
  my $i;
  my $dist;
  my ($bpxmid, $bpymid);
  for($i = 1; $i < scalar(@{$bpxbAR}); $i++) { 
    $dist = distance_between_two_points($x, $y, $bpxbAR->[$i], $bpybAR->[$i]);
    #print STDERR ("TMP x: $x y: $y dist: $dist\n");
    if($dist < $mindist) { 
      #printf STDERR ("TMP x: $x y: $y $i dist $dist $bpxbAR->[$i], $bpybAR->[$i]\n");
      $mindist = $dist; 
    }
    $dist = distance_between_two_points($x, $y, $bpxeAR->[$i], $bpyeAR->[$i]);
    if($dist < $mindist) { 
      #printf STDERR ("TMP x: $x y: $y $i dist $dist $bpxeAR->[$i], $bpyeAR->[$i]\n");
      $mindist = $dist; 
    }
    $bpxmid = ($bpxbAR->[$i] + $bpxeAR->[$i]) / 2.;
    $bpymid = ($bpybAR->[$i] + $bpyeAR->[$i]) / 2.;
    $dist = distance_between_two_points($x, $y, $bpxmid, $bpymid);
    if($dist < $mindist) { 
      #printf("TMP x: $x y: $y $i dist $dist $bpxeAR->[$i], $bpyeAR->[$i]\n");
      $mindist = $dist; 
    }
  }
  return $mindist;
}
############################################################
sub angle_between_two_points 
{
  my $deg;
  my ($x1, $y1, $x2, $y2) = @_;
  if(abs($x1-$x2) <= 0.0001) { $deg = 90.; }
  else { $deg = rad2deg(atan((abs($y1-$y2))/(abs($x1-$x2)))); }
  return $deg;
}
############################################################
sub distance_between_two_points 
{
  my $deg;
  my ($x1, $y1, $x2, $y2) = @_;
  my $dist = sqrt(((abs($x1-$x2))**2.) + ((abs($y1-$y2))**2.));
  return $dist;
}
############################################################
sub dist2closest_other_residue
{
  my ($xAR, $yAR, $posn, $x, $y) = @_;
  my $mindist = 100000;
  my $cres = 0;
  my ($i, $dist);
  for($i = 1; $i < scalar(@{$xAR}); $i++) { 
    if($i != $posn) { 
      $dist = distance_between_two_points($x, $y, $xAR->[$i]+($fontsize/2.), $yAR->[$i]+($fontsize/2.));
      if($dist < $mindist) { $mindist = $dist; $cres = $i; }
    }
  }
  return $dist;
}
##########################################################
sub input_svg_file
{
  my ($svgIFH, $clen, $nt_nAR, $nt_xAR, $nt_yAR) = @_;

  # unsorted nucleotide data
  my %nt_nH = (); # actual nucleotide
  my %nt_xH = (); # x coord
  my %nt_yH = (); # y coord
  # sorted nucleotide data
  # Input SVG file
  my $i = 0;
  my ($line, $scale, $nt, $x, $y, $id);
  while($line = <$svgIFH>) { 
    chomp $line;
    if($line =~ m/^<svg/) { 
      $scale = svg_input_scale($svg_or_crwIFH, $line)
        }
    if($line =~ m/^<text/) { 
      if(! defined $scale) { die "ERROR did not find SVG header line"; }
      svg_input_nucleotide($svg_or_crwIFH, \$nt, \$x, \$y, \$id);
      $nt_nH{$id} = $nt;
      $nt_xH{$id} = $x;
      $nt_yH{$id} = $y;
      # print("id: $id n: $nt x: $x y: $y\n");
    }
  }

  # put nucleotides in correct order
  # example order:
  # $id value        "1", "2", "3", "3.1", "4", "5", "5.1", "5.2", "6"
  # order in @nt_*A:  0    1    2    3      4    5    6      7      8
  my $nnt = scalar(keys(%nt_nH));
  my %usedH = ();
  my ($j, $ij);
  my $max_id = 0;
  foreach $id (keys %nt_nH) { 
    $usedH{$id} = 0; 
    if($id > $max_id) { $max_id = $id; }
  }
  if($max_id =~ m/\./) { # final nucleotide can't have a '.' in it, this simplifies loop below that puts them in proper order
    die "ERROR final nucleotide's id is $max_id: it must be an integer"; 
  }
  
  $i = 1;
  while($i <= $max_id) { 
    if(exists($nt_nH{$i})) { 
      push(@{$nt_nAR}, $nt_nH{$i});
      #printf("pushing i: $i (total: %d)\n", scalar(@{$nt_nAR}));
      push(@{$nt_xAR}, $nt_xH{$i} / $scale);
      push(@{$nt_yAR}, (792. - $nt_yH{$i}) / $scale);
      $usedH{$i} = 1;
      $j = 1;
      $ij = $i . "." . $j;
      while(exists($nt_nH{$ij})) { 
        push(@{$nt_nAR}, $nt_nH{$ij});
        #printf("pushing ij: $ij (total: %d)\n", scalar(@{$nt_nAR}));
        push(@{$nt_xAR}, $nt_xH{$ij} / $scale);
        push(@{$nt_yAR}, (792. - $nt_yH{$ij}) / $scale);
        $usedH{$ij} = 1;
        $j++;
        $ij = $i . "." . $j;
      }
    }
    $i++;
  }
  # sanity check: make sure we got all of them
  my $nmissed = 0;
  foreach $id (keys %nt_nH) { 
    if(! $usedH{$id}) { 
      print "problem sorting nucleotides, $id not included\n"; 
      $nmissed++;
    }
  }
  if($nmissed > 0) { 
    die "ERROR $nmissed nucleotides not sorted properly"; 
  }
  my $obs_clen = scalar(@nt_nA);
  if($obs_clen != $clen) { die "ERROR input alignment had $clen alignment positions, but read $obs_clen residues in SVG file"; }
  return;
}

############################################
sub input_crw_ps_file { 
  my ($psIFH, $clen, $nt_nAR, $nt_xAR, $nt_yAR) = @_;

  my $nnt = 0;
  my $xoff = 0;
  my $yoff = 0;
  while($line = <$psIFH>) { 
    if($line =~ /^(\-?\d+\.\d\d)\s(\-?\d+\.\d\d) translate/) { 
      $xoff += $1;
      $yoff += $2;
    }
    if($line =~ /^(\d+\.\d\d)\s(\d+\.\d\d) scale/) { 
      if($1 ne $2) { die "ERROR scale line not uniform!\n"; }
      $xoff /= $1;
      $yoff /= $2;
    }
    if($line =~ /^(\(\w\))\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)\s+lwstring/) { 
      my ($nt, $x, $y) = ($1, $2, $3);
      if($nt =~ m/[ACGU]/){ 
        push(@nt_nA, $nt);
        push(@nt_xA, $x+$xoff);
        push(@nt_yA, $y+$yoff);
        $nnt++;
      }
    }
  }
  if($nnt != $clen) { die "ERROR SS_cons length not equal to number of nts read from CRW ps file ($clen != $nnt)"; }
}
