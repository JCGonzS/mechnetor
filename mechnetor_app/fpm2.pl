#!/usr/bin/perl -w

$|=1;

# Find motif in protein sequence database

$mode = 4;
$add = 0;
$case = 0;
$dna_trans = 0;

# These are human uniref 90 frequencies = change if desired
$AA_FREQ{"A"}=0.07092642; $AA_FREQ{"B"}=0.00001075; $AA_FREQ{"C"}=0.02441919; $AA_FREQ{"D"}=0.04412000;
$AA_FREQ{"E"}=0.06601801; $AA_FREQ{"F"}=0.03535145; $AA_FREQ{"G"}=0.07044761; $AA_FREQ{"H"}=0.02679318;
$AA_FREQ{"I"}=0.04067031; $AA_FREQ{"K"}=0.05324042; $AA_FREQ{"L"}=0.09892366; $AA_FREQ{"M"}=0.02183118;
$AA_FREQ{"N"}=0.03288463; $AA_FREQ{"P"}=0.06678249; $AA_FREQ{"Q"}=0.04686112; $AA_FREQ{"R"}=0.05886737;
$AA_FREQ{"S"}=0.08569498; $AA_FREQ{"T"}=0.05519229; $AA_FREQ{"U"}=0.00000109; $AA_FREQ{"V"}=0.06000912;
$AA_FREQ{"W"}=0.01485423; $AA_FREQ{"X"}=0.00090436; $AA_FREQ{"Y"}=0.02518633; $AA_FREQ{"Z"}=0.00000981;

$codon_file = "dbs/codon.txt";

for($i=0; $i<=$#ARGV; ++$i) {
   if($ARGV[$i] eq "-m") { # Mode
        if(($i+1)>($#ARGV)) { exit_error(); }
        $mode = $ARGV[$i+1];
        $i++;
   } elsif($ARGV[$i] eq "-add") { # no gaps
        if(($i+1)>($#ARGV)) { exit_error(); }
        $add  = $ARGV[$i+1];
        $i++;
   } elsif($ARGV[$i] eq "-case") { #
        $case=1;
   } elsif($ARGV[$i] eq "-dna_trans") { #
        $dna_trans=1;
   } elsif(substr($ARGV[$i],0,1) eq "-") {
        print STDERR "Option $ARGV[$i] not recognised\n\n";
        exit_error();
   } else {
       $motif = $ARGV[$i];
       if(length($motif)<1) {
        print STDERR "Motif is too short or missing\n\n";
        exit_error();
       }
   }
}

if($dna_trans==1) {
 open(IN,$codon_file) || die "Error reading $codon_file\n";
 while(<IN>) {
  chomp;
  @T = split(/\s+/);
  for($i=1; $i<=$#T; ++$i) {
   $CODON2AA{$T[$i]} = $T[0];
#   $AA2CODON{$T[0]}{$T[$i]}=1;
   $T[$i] =~ s/U/T/g;
   $CODON2AA{$T[$i]} = $T[0];
#   $AA2CODON{$T[0]}{$T[$i]}=1;
  }
 }
 close(IN);
}


$pmotif = 1;
for($i=0; $i<length($motif); ++$i) {
  $aa = substr($motif,$i,1);
  if(defined($AA_FREQ{$aa})) { $pmotif *= $AA_FREQ{$aa} }
}
if(($pmotif < 0.001) && ($pmotif>0)) {
    $smotif = sprintf("%9.3e",$pmotif);
} else {
     $smotif = sprintf("%9.6f",$pmotif);
}
$dmotif = $motif;
$motif =~ s/x/./g;

#printf(STDERR "Background probability calculated as %10.8e\n",$pmotif);

sub exit_error {
   die "fpm.pl [motif] < [fasta sequence file]\n
  Options:
   -m [0,1,2]     -> (mode) 0: annotate sequences with instances (default)
                  ->        1: show only instances
                  ->        2: silent (just output summary at end)
                  ->        3: show only instances (table form)
                  ->        4: as 0 but only showing those containing sequences
   -case          -> show sequence outside of the motif as lower case
                     (mode 0 only ; by default is off)
   -add [integer] -> number of additional residues to add to N- C- term
                     (mode 1 only ; default is 0)\n";
}

if(!defined($motif)) {
  print STDERR "Motif is too short or missing\n\n";
  exit_error();
}
$n_seq=0;
$n_occ=0;
$t_seq=0;
$t_aas=0;
while(<stdin>) {
   if(/>/) {
      if($t_seq>0) {
        compare_seq_multiple();  # Do the last sequence
      }
      $t_seq++;
      chomp;
      $this_entry_id = $_;
      if($this_entry_id !~ / /) { $this_entry_id .= " " }
      $this_entry_seq = "";
      $seq = "";
      $gn = "no-gn-found";
      if(/GN=/) {
       $gn = $_;
       chomp $gn;
       $gn =~ s/.*GN=//;
       $gn =~ s/ .*//;
      }
#      $GN{$this_entry_id}=1;
   } else {
      $this_entry_seq .= $_;
      chomp;
      $seq .= uc($_);
   }

}
if($t_seq>0) {
   compare_seq_multiple(); # Do the last sequence
   $f_seq = $n_seq/$t_seq;
   if(($f_seq < 0.001) && ($f_seq>0)) {
        $s1 = sprintf("%9.3e",$f_seq);
   } else {
        $s1 = sprintf("%9.6f",$f_seq);
   }
   $f_occ = $n_occ/$t_aas;
   if(($f_occ < 0.001) && ($f_occ>0)) {
        $s2 = sprintf("%9.3e",$f_occ);
   } else {
        $s2 = sprintf("%9.6f",$f_occ);
   }
  #  printf(STDERR "M: %10s SEQ: %4d / %5d = %s AAS: %4d / %8d = %s Exp AA prob: %s Ratio: %8.4f\n",
  #      $dmotif,
  #      $n_seq,$t_seq,$s1,
  #      $n_occ,$t_aas,$s2,$smotif,$f_occ/$pmotif);
}


sub compare_seq_multiple {
   $temp_seq = $seq; $temp_seq =~ s/[Xx]//g;
   $t_aas += length($temp_seq);
   if($seq =~ /$motif/) { $n_seq++; }
   $output_s = "";
   $tmp_seq = $seq;
   if((($mode==0) || ($mode==4)) && ($case==1)) { $seq = lc($seq); }
   $n_instances=0;
   while(uc($tmp_seq) =~ /$motif/) {
#    print "HERE $n_seq ",length($tmp_seq),"\n",substr($tmp_seq,0,1000),"\n";
     $n_instances++;
     $n_occ++;
     $s = length($`)+1;
     $e = length($seq)-length($');
     $l = $e - $s+1;
     if($l<0) { print STDERR "ERROR funning start/end\n"; last; }
     $output_s .= sprintf(" %4d - %4d; ",$s,$e);
     $suffix = sprintf("/%d-%d",$s,$e);
     $tmp_seq2 = substr($tmp_seq,0,$s-1);
     for($i=0; $i<$l; ++$i) { $tmp_seq2 .= "x"; }
     $tmp_seq2 .= substr($tmp_seq,$e);
#     print "$seq\n$tmp_seq2\n";
     $tmp_seq = $tmp_seq2;
#     printf("Found motif at %4d - %4d (%2d) %s\n",$s,$e,$l,substr($tmp_seq,0,1000));
      if((($mode==0) || ($mode==4)) && ($case==1)) {
         $tmp_seq2 = substr($seq,0,$s-1);
         $tmp_seq2 .= uc(substr($seq,$s-1,$l));
         $tmp_seq2 .= substr($seq,$e);
         if(length($seq) != length($tmp_seq2)) { die "error in length\n"; }
         $seq = $tmp_seq2;
      } elsif(($mode==1) || ($mode==3)) {
         $s-= $add; if($s<0) { $s=0; }
         $e+= $add; if($e>length($seq)) { $e=length($seq); }
         $l = $e - $s +1;
         $ts = $this_entry_id; $ts =~ s/ /$suffix /;
         if($mode==1) {
            printf("%-25s %4d - %4d (%2d)\n%s\n",
               $ts,$s,$e,$l,substr($seq,$s-1,$l));
         } else {
            $ts =~ s/ .*//;
            printf("%-25s %10s %4d - %4d (%2d) : %s",
               $ts,$gn,$s,$e,$l,substr($seq,$s-1,$l));
            if($dna_trans==1) {
              # Assume cDNA so just report any overlapping AAs
              $s2 = $s - 2;
              $e2 = $e;
              if($s2 < 1) { $s2 = 1 }
              if($e2 > length($seq)) { $e2 = length($seq) }
              for($k=$s2; $k<=$e2; ++$k) {
               if((($k-1)%3)==0) {
                $codon = substr($seq,$k-1,3);
                if(length($codon)==3) {
                 if(defined($CODON2AA{$codon})) {
                  $trans_AA = $CODON2AA{$codon};
                  $trans_AA_pos = int(($k-1)/3)+1;
                  if(!defined($ALREADY{$trans_AA}{$trans_AA_pos})) { printf(" %s%d",$trans_AA,$trans_AA_pos); }
                  $ALREADY{$trans_AA}{$trans_AA_pos}=1;
                 } else {
                  print "No AA defined for codon $codon\n";
                 }
                }
               }
              }
            }
            print "\n";
         }
      }
   }
   if(($mode==3) && ($n_instances==0)) {
      $ts = $this_entry_id;  $ts =~ s/ .*//;
      printf("%-37s %10s : no hits\n",$ts,$gn);
   }
   if(($mode==0) || (($mode==4) && ($n_instances>0))) {
       print $this_entry_id," ;  Motif: ",$dmotif," is found ",$n_instances," times ",$output_s,"\n";
#       print $this_entry_seq;
        for($i=0; $i<length($seq); ++$i) {
          print substr($seq,$i,1);
          if((($i+1)%60)==0) { print "\n"; }
        }
        print "\n"

   }
}
