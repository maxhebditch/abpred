#!/usr/bin/perl -w

$file_here		= 'blah.txt';
open (BLAH, ">>$file_here");		# file used for comments throughout server code

$window[1]	= 21;
$window_half[1]	= 10;
$window[2]	= 51;
$window_half[2]	= 25;
$window_flag[1]	= "WIN21";
$window_flag[2]	= "WIN51";

@aa_single	= ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");

if ((exists $ENV{MEMBRANE_EXCLUDE}) and ($ENV{MEMBRANE_EXCLUDE} eq "YES")) {
  $mem_exclude	= "yes";
  $KD_threshold	= 1.6;	# value of 19 win avg (-4.5 to 4.5) KD suggested by KD - we'll use for 21 win
  printf BLAH "membrane proteins will be excluded, based on KD of 21 aa window > ?? value absent \n";
} else {
  $mem_exclude	= "no";
  $KD_threshold	= 1.6;	# value of 19 win avg (-4.5 to 4.5) KD suggested by KD - we'll use for 21 win, still use even if not excl
}

# e.g. for excel correlation calc, if we have the abundance/solub data in the ID line
if ((exists $ENV{WHOLE_SEQ_OUT}) and ($ENV{WHOLE_SEQ_OUT} eq "yes")) {
  $whole_seq_out	= "yes";
  printf BLAH "WHOLE_SEQ_OUT set: for WHOLE-SEQ only, output file ID line with whole seq data, in whole_seq_data.txt\n";
} else {
  $whole_seq_out	= "no";
}


$file_here		= 'composition.in';
open (COMP_IN, "<$file_here") or die "cannot open composition.in\n";
$nseq		= 0;
$first_seq		= 'yes';
while ($line 		= <COMP_IN>) {
  chomp $line;
  @words 	= split ("",$line);
  @fields	= split (" ",$line);
  if (exists $words[0]) {
    if ($words[0] eq ">") {				# new protein ID
      if ($first_seq eq 'yes') {			# first ID, no seqbuf
        $seqbuf		= "";
	$idbuf		= $fields[0];
	$textbuf	= $line;
  	$first_seq	= 'no';
      } else {
        if ($seqbuf ne "") {				# only process further if non-null seq
	  $nseq++;
	  $seqtext[$nseq]		= $textbuf;
	  $seqid[$nseq]			= $idbuf;
	  $seqtmp{$seqid[$nseq]}	= $seqbuf;	# store away the sequence
	}
	$idbuf		= $fields[0];
	$textbuf	= $line;
	$seqbuf		= "";
      }
    } else {
      if ($first_seq eq 'yes') {
	printf BLAH "did not see first seq ID\n";
  	exit;
      }
      $seqbuf	.= $line;
    }
  }
}

if ($seqbuf ne "") {				# catch the last seq, if non-null
  $nseq++;
  $seqtext[$nseq]		= $textbuf;
  $seqid[$nseq]			= $idbuf;
  $seqtmp{$seqid[$nseq]}	= $seqbuf;	# store away the sequence
}

$nseqs			= $nseq;
close (COMP_IN);

$file_here		= 'composition_all.out';
open (COMP_OUT, ">$file_here") or die "cannot open composition_all.out\n";
if ($whole_seq_out eq "yes") {
  $file_here		= 'whole_seq_data.txt';
  open (DATA_OUT, ">$file_here") or die "cannot open whole_seq_data.txt\n";
}

printf COMP_OUT "\nALL WIN21-MAX, WIN21-MIN, WIN51-MAX, WIN51-MIN data lines for each sequence\n";
printf COMP_OUT "output header lines flagged by each data type for grep purposes\n\n";

@flags = ('WHOLE-SEQ', 'WIN21-MAX', 'WIN21-MIN', 'WIN51-MAX', 'WIN51-MIN');
foreach $win_flag (@flags) {
  printf COMP_OUT "$win_flag, ORF-ID,K-R,D-E,naa,totperc,";
  printf COMP_OUT "A,C,D,E,F,G,H,I,K,L,";
  printf COMP_OUT "M,N,P,Q,R,S,T,V,W,Y,";
  printf COMP_OUT "K+R,D+E,K+R-D-E,K+R+D+E,F+W+Y,pI,";
  printf COMP_OUT "KyteDoo,abs-charge,FoldIndex,disorder,entropy,betapropensity\n\n";
  if (($whole_seq_out eq 'yes') and ($win_flag eq 'WHOLE-SEQ')) {
    printf DATA_OUT "$win_flag, ORF-ID,K-R,D-E,naa,totperc,";
    printf DATA_OUT "A,C,D,E,F,G,H,I,K,L,";
    printf DATA_OUT "M,N,P,Q,R,S,T,V,W,Y,";
    printf DATA_OUT "K+R,D+E,K+R-D-E,K+R+D+E,F+W+Y,pI,";
    printf DATA_OUT "KyteDoo,abs-charge,FoldIndex,disorder,entropy,betapropensity\n";
  }
}

$Kall		= 0;
$Rall		= 0;
$Dall		= 0;
$Eall		= 0;
$Fall		= 0;
$Wall		= 0;
$Yall		= 0;
$naa_all	= 0;
$nmem_exclude	= 0;
$nall_exclude	= 0;
for ($nseq=1; $nseq<=$nseqs; $nseq++) {
	$seqhere	= $seqtmp{$seqid[$nseq]};
	$lentmp		= length $seqhere;

	$KD_max_original = 999;

	@words		= split ("",$seqhere);
	$naa		= 0;
	$Kcount		= 0;
	$Rcount		= 0;
	$Dcount		= 0;
	$Ecount		= 0;
	$Fcount		= 0;
	$Wcount		= 0;
	$Ycount		= 0;
	for ($aa=1; $aa<=20; $aa++) {
	  $aa_count[$aa-1]	= 0;
	}
	@seq_for_pI	= ();			# hopefully this zeroes the array
	$nseq_for_pI	= 0;
	for ($n=1; $n<=$lentmp; $n++) {
	  if ($words[$n-1] ne " ") {
	    if ($words[$n-1] eq "*") {
	      if ($naa == 0) {
	        printf BLAH "STOPPING - first seq char is a * - Plation\n";
		exit;
	      }
	      $phosupd[$naa]	= "P";
	    } else {
	      $seq_for_pI[$nseq_for_pI]	= $words[$n-1];
	      $nseq_for_pI++;
	      $naa++;
	      $sequpd[$naa]	= $words[$n-1];
	      $phosupd[$naa]	= " ";
	      if ($sequpd[$naa] eq 'K') { $Kcount++; }
	      if ($sequpd[$naa] eq 'R') { $Rcount++; }
	      if ($sequpd[$naa] eq 'D') { $Dcount++; }
	      if ($sequpd[$naa] eq 'E') { $Ecount++; }
	      if ($sequpd[$naa] eq 'F') { $Fcount++; }
	      if ($sequpd[$naa] eq 'W') { $Wcount++; }
	      if ($sequpd[$naa] eq 'Y') { $Ycount++; }
	      for ($aa=1; $aa<=20; $aa++) {
	        if ($sequpd[$naa] eq $aa_single[$aa-1]) { $aa_count[$aa-1]++; }
	      }
	    }
	  }
	}

        $pI_all		= protein_pI (@seq_for_pI);
	$null_return	= 'no';
	$KD_all		= protein_KD (@seq_for_pI);
	if ($KD_all == 999) { $null_return = 'yes'; }
	$Q_all		= protein_Q (@seq_for_pI);
	if ($Q_all == 999) { $null_return = 'yes'; }
	$absQ_all	= abs ($Q_all);
        $foldind_all	= 2.785*$KD_all - $absQ_all - 1.151;
	$disorder_all	= protein_disorder (@seq_for_pI);
	if ($disorder_all == 999) { $null_return = 'yes'; }
	$ent_all	= protein_ent (@seq_for_pI);
	if ($ent_all == 999) { $null_return = 'yes'; }
	$betastrand_all	= protein_betastrand (@seq_for_pI);
	if ($betastrand_all == 999) { $null_return = 'yes'; }

	$perc_total		= 0;		# work out and output the various quantities for whole sequence
	for ($aa=1; $aa<=20; $aa++) {
	  $aa_perc[$aa-1]	= (100*$aa_count[$aa-1]/$naa);
	  $perc_total		+= $aa_perc[$aa-1];
	}
	if ($naa != 0) {
	  $Kcomp	= $Kcount/$naa;
	  $Rcomp	= $Rcount/$naa;
	  $Dcomp	= $Dcount/$naa;
	  $Ecomp	= $Ecount/$naa;
	  $Fcomp	= $Fcount/$naa;
	  $Wcomp	= $Wcount/$naa;
	  $Ycomp	= $Ycount/$naa;
	  $KmR_perc	= 100*($Kcomp - $Rcomp);
	  $DmE_perc	= 100*($Dcomp - $Ecomp);
	} else {
	  $Kcomp	= 0;			# not really 0 of course
	  $Rcomp	= 0;			# also not really 0
	  $Dcomp	= 0;			# not really 0 of course
	  $Ecomp	= 0;			# also not really 0
	  $Fcomp	= 0;
	  $Wcomp	= 0;
	  $Ycomp	= 0;
	  $KmR_perc	= 0;
	  $DmE_perc	= 0;
	}
	$KR_perc	= 100*($Kcomp+$Rcomp);
	$DE_perc	= 100*($Dcomp+$Ecomp);
	$FWY_perc	= 100*($Fcomp+$Wcomp+$Ycomp);
	$KR_DE_perc	= $KR_perc - $DE_perc;
	$KRplusDE_perc	= $KR_perc + $DE_perc;
	$Kall		+= $Kcount;
	$Rall		+= $Rcount;
	$Dall		+= $Dcount;
	$Eall		+= $Ecount;
	$Fall		+= $Fcount;
	$Wall		+= $Wcount;
	$Yall		+= $Ycount;
	$naa_all	+= $naa;

    if ($null_return ne 'no') {
	  printf COMP_OUT "\nSTOPPING - PROBLEM,$seqid[$nseq],null_return\n";
	  exit;
	}

	if ($naa != 0) {			# now for the 21 and 51 aa windows, maxima and minima
	  for ($nwin=1; $nwin<=2; $nwin++) {	# 1 = 21 and 2 = 51
	    $window2	= $window_half[$nwin];
	    $win_flag	= $window_flag[$nwin];
	    $at_least_one	= "no";		# overwritten if we see a full length window

	    $KR_perc_max 	= -100;
	    $KR_perc_min 	= 100;
	    $DE_perc_max 	= -100;
	    $DE_perc_min 	= 100;
	    $KR_DE_perc_max 	= -100;
	    $KR_DE_perc_min 	= 100;
	    $KRplusDE_perc_max 	= -100;
	    $KRplusDE_perc_min 	= 100;
	    $KmR_perc_max	= -100;
	    $KmR_perc_min	= 100;
	    $DmE_perc_max	= -100;
	    $DmE_perc_min	= 100;
	    $FWY_perc_max	= -100;
	    $FWY_perc_min	= 100;
	    for ($aa=1; $aa<=20; $aa++) {
	      $aa_perc_min[$aa-1]	= 0;
	      $aa_perc_max[$aa-1]	= 0;
	    }
	    $pI_max		= 1;
	    $pI_min		= 14;

	    $KD_max		=  0;
	    $KD_min		=  1;

	    $absQ_max		=  0;
	    $absQ_min		=  1;

	    $foldind_max	= -2;
	    $foldind_min	=  2;

	    $disorder_max	= -1;
	    $disorder_min	=  1;

	    $ent_max		=  0;
	    $ent_min		=  5;

	    $betastrand_max	= -1;
	    $betastrand_min	=  1;

	    for ($maa=1; $maa<=$naa; $maa++) {
	      $Kwin		= 0;
	      $Rwin		= 0;
	      $Dwin		= 0;
	      $Ewin		= 0;
	      $Fwin		= 0;
	      $Wwin		= 0;
	      $Ywin		= 0;
	      for ($aa=1; $aa<=20; $aa++) {
	        $aa_count_win[$aa-1]	= 0;
	      }

	      $win_start	= $maa - $window2;
	      if ($win_start < 1) {$win_start = 1; }
	      $win_end		= $maa + $window2;
	      if ($win_end > $naa) {$win_end = $naa; }
	      @seq_for_pI	= ();			# hopefully this zeroes the array
	      $nseq_for_pI	= 0;
	      for ($win=$win_start; $win<=$win_end; $win++) {
	        if ($sequpd[$win] eq 'K') { $Kwin++; }
	        if ($sequpd[$win] eq 'R') { $Rwin++; }
	        if ($sequpd[$win] eq 'D') { $Dwin++; }
	        if ($sequpd[$win] eq 'E') { $Ewin++; }
	        if ($sequpd[$win] eq 'F') { $Fwin++; }
	        if ($sequpd[$win] eq 'W') { $Wwin++; }
	        if ($sequpd[$win] eq 'Y') { $Ywin++; }
	        for ($aa=1; $aa<=20; $aa++) {
	          if ($sequpd[$win] eq $aa_single[$aa-1]) { $aa_count_win[$aa-1]++; }
	        }
	        $seq_for_pI[$nseq_for_pI]	= $sequpd[$win];
	        $nseq_for_pI++;
	      }
	      $num_win		= $win_end - $win_start + 1;
	      if ($num_win == $window[$nwin]) {			# process
	        $at_least_one		= "yes";		# we have a full length window
		$perc_total_win		= 0;
		for ($aa=1; $aa<=20; $aa++) {
	  	  $aa_perc_win[$aa-1]	= (100*$aa_count_win[$aa-1]/$num_win);
	  	  $perc_total_win	+= $aa_perc_win[$aa-1];
		}
	        $Kcomp_win	= $Kwin/$num_win;
	        $Rcomp_win	= $Rwin/$num_win;
	        $Dcomp_win	= $Dwin/$num_win;
	        $Ecomp_win	= $Ewin/$num_win;
		$Fcomp_win	= $Fwin/$num_win;
		$Wcomp_win	= $Wwin/$num_win;
		$Ycomp_win	= $Ywin/$num_win;
	        $KmR_perc_win	= 100*($Kcomp_win - $Rcomp_win);
	        $DmE_perc_win	= 100*($Dcomp_win - $Ecomp_win);
		$KR_perc_win	= 100*($Kcomp_win+$Rcomp_win);
		$DE_perc_win	= 100*($Dcomp_win+$Ecomp_win);
		$KR_DE_perc_win	= $KR_perc_win - $DE_perc_win;
		$KRplusDE_perc_win	= $KR_perc_win + $DE_perc_win;
		$FWY_perc_win	= 100*($Fcomp_win+$Wcomp_win+$Ycomp_win);

		if ($KR_perc_win > $KR_perc_max) { $KR_perc_max = $KR_perc_win; }
		if ($KR_perc_win < $KR_perc_min) { $KR_perc_min = $KR_perc_win; }
		if ($DE_perc_win > $DE_perc_max) { $DE_perc_max = $DE_perc_win; }
		if ($DE_perc_win < $DE_perc_min) { $DE_perc_min = $DE_perc_win; }
		if ($KR_DE_perc_win > $KR_DE_perc_max) { $KR_DE_perc_max = $KR_DE_perc_win; }
		if ($KR_DE_perc_win < $KR_DE_perc_min) { $KR_DE_perc_min = $KR_DE_perc_win; }
		if ($KRplusDE_perc_win > $KRplusDE_perc_max) { $KRplusDE_perc_max = $KRplusDE_perc_win; }
		if ($KRplusDE_perc_win < $KRplusDE_perc_min) { $KRplusDE_perc_min = $KRplusDE_perc_win; }
		if ($KmR_perc_win > $KmR_perc_max) { $KmR_perc_max = $KmR_perc_win; }
		if ($KmR_perc_win < $KmR_perc_min) { $KmR_perc_min = $KmR_perc_win; }
		if ($DmE_perc_win > $DmE_perc_max) { $DmE_perc_max = $DmE_perc_win; }
		if ($DmE_perc_win < $DmE_perc_min) { $DmE_perc_min = $DmE_perc_win; }
		if ($FWY_perc_win > $FWY_perc_max) { $FWY_perc_max = $FWY_perc_win; }
		if ($FWY_perc_win < $FWY_perc_min) { $FWY_perc_min = $FWY_perc_win; }

	        for ($aa=1; $aa<=20; $aa++) {
	          if ($aa_perc_win[$aa-1] > $aa_perc_max[$aa-1]) { $aa_perc_max[$aa-1] = $aa_perc_win[$aa-1]; }
	          if ($aa_perc_win[$aa-1] < $aa_perc_min[$aa-1]) { $aa_perc_min[$aa-1] = $aa_perc_win[$aa-1]; }
	        }

                $pI_win		= protein_pI (@seq_for_pI);
		if ($pI_win > $pI_max) { $pI_max = $pI_win; }
		if ($pI_win < $pI_min) { $pI_min = $pI_win; }

		$null_return	= 'no';
		$KD_win		= protein_KD (@seq_for_pI);
		if ($KD_win == 999) { $null_return = 'yes'; }
		$Q_win		= protein_Q (@seq_for_pI);
		if ($Q_win == 999) { $null_return = 'yes'; }
		$absQ_win	= abs ($Q_win);
	        $foldind_win	= 2.785*$KD_win - $absQ_win - 1.151;
		$disorder_win	= protein_disorder (@seq_for_pI);
		if ($disorder_win == 999) { $null_return = 'yes'; }
		$ent_win	= protein_ent (@seq_for_pI);
		if ($ent_win == 999) { $null_return = 'yes'; }
		$betastrand_win	= protein_betastrand (@seq_for_pI);
		if ($betastrand_win == 999) { $null_return = 'yes'; }

		if ($KD_win > $KD_max) {
		  $KD_max = $KD_win;
		  $KD_max_original	= (9*$KD_max) - 4.5	# rescale back to KD proper, use if we're exlcuding mem proteins
		}
		if ($KD_win < $KD_min) { $KD_min = $KD_win; }
		if ($absQ_win > $absQ_max) { $absQ_max = $absQ_win; }
		if ($absQ_win < $absQ_min) { $absQ_min = $absQ_win; }
		if ($foldind_win > $foldind_max) { $foldind_max = $foldind_win; }
		if ($foldind_win < $foldind_min) { $foldind_min = $foldind_win; }
		if ($disorder_win > $disorder_max) { $disorder_max = $disorder_win; }
		if ($disorder_win < $disorder_min) { $disorder_min = $disorder_win; }
		if ($ent_win > $ent_max) { $ent_max = $ent_win; }
		if ($ent_win < $ent_min) { $ent_min = $ent_win; }
		if ($betastrand_win > $betastrand_max) { $betastrand_max = $betastrand_win; }
		if ($betastrand_win < $betastrand_min) { $betastrand_min = $betastrand_win; }

	      }
	    }		# end of this (21 51) window sliding along for this sequence

	    if ($nwin == 1) {			# 21 aa window
#	      if ($mem_exclude eq "yes") {	# put lower down to give flag pred TM, even if not skipping it
	        if ($win_flag ne "WIN21") {
		  printf BLAH "STOPPING - win_flag should be 21 for KD-membrane protein check\n";
		  exit;
		}
	        if ($null_return ne 'no') { goto SKIPPED; }	# skip protein since we cannot KD test a null return
		if ($at_least_one ne "yes") { goto SKIPPED; }	# skip protein since no window to KD test
		if ($KD_max_original == 999) {
		  printf BLAH "STOPPING - KD_max_original unset\n";
		} else {
		  if ($KD_max_original > $KD_threshold) {	# any KD looks like mem window, report, count and outta here
	            if ($mem_exclude eq "yes") {
		      $nmem_exclude++;
		      printf BLAH "MEM EXCLUDE - $seqid[$nseq] - with KD, threshold = $KD_max_original $KD_threshold\n";
		      goto SKIPPED;
		    } else {
		      printf BLAH "\nTM REGION PREDICTED for $seqid[$nseq] - with non-normalised KD of ";
		      printf BLAH "%6.2f", $KD_max_original;
		      printf BLAH " compared with threshold of $KD_threshold\n";
		    }
		  }
		}
	      printf COMP_OUT "\n";
	      printf COMP_OUT "WHOLE-SEQ,$seqid[$nseq],";
	      printf COMP_OUT "%7.2f,%7.2f,%5.0f,%7.2f", $KmR_perc, $DmE_perc, $naa, $perc_total;
	      for ($aa=1; $aa<=20; $aa++) { printf COMP_OUT ",%6.2f", $aa_perc[$aa-1]; }
	      printf COMP_OUT ",%7.2f,%7.2f,%7.2f,%7.2f,%7.2f,%6.2f", $KR_perc, $DE_perc, $KR_DE_perc, $KRplusDE_perc, $FWY_perc, $pI_all;
	      printf COMP_OUT ",%8.3f,%8.3f,%8.3f,%8.3f, %8.3f, %8.3f", $KD_all, $absQ_all, $foldind_all, $disorder_all, $ent_all, $betastrand_all;
	      printf COMP_OUT "\n";

	      if ($whole_seq_out eq 'yes') {
#	        printf DATA_OUT "\n";
	        printf DATA_OUT "WHOLE-SEQ,$seqtext[$nseq],";		# seqtext (rather than seqid) should be the whole ID line
	        printf DATA_OUT "%7.2f,%7.2f,%5.0f,%7.2f", $KmR_perc, $DmE_perc, $naa, $perc_total;
	        for ($aa=1; $aa<=20; $aa++) { printf DATA_OUT ",%6.2f", $aa_perc[$aa-1]; }
	        printf DATA_OUT ",%7.2f,%7.2f,%7.2f,%7.2f,%7.2f,%6.2f", $KR_perc, $DE_perc, $KR_DE_perc, $KRplusDE_perc, $FWY_perc, $pI_all;
	        printf DATA_OUT ",%8.3f,%8.3f,%8.3f,%8.3f, %8.3f, %8.3f", $KD_all, $absQ_all, $foldind_all, $disorder_all, $ent_all, $betastrand_all;
	        printf DATA_OUT "\n";
	      }

	    }		# end nwin=1 check

	    if ($null_return eq 'no') {		# at present have no reporting of null returns for windows
	     if ($at_least_one eq "yes") {	# we have at least one full length window
						# OTHERWISE we simply MISS OUT this write
						# so NOT ALL SEQS will have writes for both windows
	      $flag_here	= $win_flag."-MAX";
	      printf COMP_OUT "$flag_here,$seqid[$nseq],";
	      printf COMP_OUT "%7.2f,%7.2f,%5.0f,%7.2f", $KmR_perc_max, $DmE_perc_max, $window[$nwin], $perc_total_win;
	      for ($aa=1; $aa<=20; $aa++) { printf COMP_OUT ",%6.2f", $aa_perc_max[$aa-1]; }
	      printf COMP_OUT ",%7.2f,%7.2f,%7.2f,%7.2f,%7.2f,%6.2f", $KR_perc_max, $DE_perc_max, $KR_DE_perc_max, $KRplusDE_perc_max, $FWY_perc_max, $pI_max;
	      printf COMP_OUT ",%8.3f,%8.3f,%8.3f,%8.3f, %8.3f, %8.3f", $KD_max, $absQ_max, $foldind_max, $disorder_max, $ent_max, $betastrand_max;
	      printf COMP_OUT "\n";
	      $flag_here	= $win_flag."-MIN";
	      printf COMP_OUT "$flag_here,$seqid[$nseq],";
	      printf COMP_OUT "%7.2f,%7.2f,%5.0f,%7.2f", $KmR_perc_min, $DmE_perc_min, $window[$nwin], $perc_total_win;
	      for ($aa=1; $aa<=20; $aa++) { printf COMP_OUT ",%6.2f", $aa_perc_min[$aa-1]; }
	      printf COMP_OUT ",%7.2f,%7.2f,%7.2f,%7.2f,%7.2f,%6.2f", $KR_perc_min, $DE_perc_min, $KR_DE_perc_min, $KRplusDE_perc_min, $FWY_perc_min, $pI_min;
	      printf COMP_OUT ",%8.3f,%8.3f,%8.3f,%8.3f, %8.3f, %8.3f", $KD_min, $absQ_min, $foldind_min, $disorder_min, $ent_min, $betastrand_min;
	      printf COMP_OUT "\n";
	     }
	    }
	  }	# end nwin = 1, 2 loop
SKIPPED:
	  $nall_exclude++;
	}
}
$KmR_all_perc	= 100*(($Kall - $Rall)/$naa_all);	# a few things over the whole sequence set
$DmE_all_perc	= 100*(($Dall - $Eall)/$naa_all);
$Kall_perc	= 100*($Kall/$naa_all);
$Rall_perc	= 100*($Rall/$naa_all);
$Dall_perc	= 100*($Dall/$naa_all);
$Eall_perc	= 100*($Eall/$naa_all);
# $FWYall_perc	= 100*(($Fall + $Wall + $Yall)/$naa_all);
printf COMP_OUT "\n";
printf COMP_OUT "SUMS, for this set of seqs, overall K minus R and D minus E percentages = ";
printf COMP_OUT "%6.3f %6.3f\n\n", $KmR_all_perc, $DmE_all_perc;
printf COMP_OUT "SUMS, for this set of seqs, overall K R D E percentages = ";
printf COMP_OUT "%6.4f %6.4f %6.4f %6.4f\n", $Kall_perc, $Rall_perc, $Dall_perc, $Eall_perc;
if ($mem_exclude eq "yes") {
  printf COMP_OUT "\n\nMEMBRANE_EXCLUDE YES set with nmem_excluded = $nmem_exclude and nall_excluded = $nall_exclude\n";
  printf BLAH "\n\nMEMBRANE_EXCLUDE yes set with nmem_excluded = $nmem_exclude and nall_excluded = $nall_exclude\n";
  if ($whole_seq_out eq 'yes') {
    printf DATA_OUT "\n\nMEMBRANE_EXCLUDE YES set with nmem_excluded = $nmem_exclude and nall_excluded = $nall_exclude\n";
  }
}

close (COMP_OUT);
if ($whole_seq_out eq 'yes') {
  close (DATA_OUT);
}

close (BLAH);

exit;


#-------------------------------------------------------------------------
#----------------------
# SUBROUTINE protein_pI
#----------------------
# pass a list of single letter amino acid codes, return the pI (pH at
# neutral charge) for this protein

sub protein_pI {

	my @aa_seq;
	my $aa;
	my %aapK;
	my %aafactor;
	my %aanum;
	my @aa_keys;
	my $pI;
	my $protein_charge;
	my $protein_charge_last;
	my $fac;
	my $pH;
	my $pH_dec;
	my $loop;
	my $pH_pK;
	my $charge_inc;

# put the passed sequence into the aa_seq array

	@aa_seq = @_;

# data initialisation --------------------------------------------------
# set up the ionisable group pKs, note that we are leaving out Cys
# (assume that cysteines are oxidised in disulphides)
# phosphates (X) are -1 always and -1 to -2 with pK 7

	%aapK = (
 	"D" =>  4.0,
	"E" =>  4.4,
 	"H" =>  6.3,
 	"K" => 10.4,
 	"R" => 12.0,
	"X" =>  7.0,
	"Y" => 10.1);

# set up the -1 or +1 factor that gives us the charge sign (acid=-,
# base=+), and is used to get the right sign of (pH-pK)

	%aafactor = (
 	"D" =>  -1.0,
	"E" =>  -1.0,
 	"H" =>   1.0,
 	"K" =>   1.0,
 	"R" =>   1.0,
	"X" =>  -1.0,
	"Y" =>  -1.0);

# initialise counter of amino acid types (for those with pKs of interest)

	%aanum = (
 	"D" =>  0,
	"E" =>  0,
 	"H" =>  0,
 	"K" =>  0,
 	"R" =>  0,
	"X" =>  0,
	"Y" =>  0);
# -----------------------------------------------------------------------

# add the number of each ionisable aa type of interest

	my $ionisable_aas = 0;
	foreach $aa (@aa_seq) {
		if (exists $aapK{$aa}) {
			$aanum{$aa}++;
			$ionisable_aas++;
		}
	}
	@aa_keys = keys %aanum;

# for any ionisable group, using the aafactor of -1 or +1 as set above,
# the charge at pH is:	(aafactor) x { 1/(1+10**[(aafactor)(pH-pK)]) }
# and the overall protein charge is the sum of such over ionisable groups

# our strategy is to come down from a high pH, look for charge sign
# swap (- to +), and then linearly interpolate, (as in the 2CPS notes)

	$loop = 0;
	$pH = 14;
	$pH_dec = 0.5;
	$protein_charge = -100;		# dummy value, calculate as we go
	while (($protein_charge < 0) and ($pH > 0.9)) {
		$loop++;
		$protein_charge_last = $protein_charge;
		$protein_charge = 0.0;
		foreach $aa (@aa_keys) {
			if ($aanum{$aa} > 0) {
				$fac = $aafactor{$aa};
				$pH_pK = $pH - $aapK{$aa};
				if ($aa eq "X") {
					$charge_inc = -1.0 + $fac * (1/(1+10**($fac*$pH_pK)));
				} else {
					$charge_inc = $fac * (1/(1+10**($fac*$pH_pK)));
				}
				$protein_charge += $aanum{$aa}*$charge_inc;
			}
		}
		$pH -= $pH_dec;	# i.e. $pH - $pH-$pH_dec
	}

# There are various conditions to look out for and signal on return:
# 	with no ionisable groups, neutral at all pH:	return 999 - NO, just return pH 7.0
#	if pI > 14 (e.g. only basic aas):		return 14
#	if pI < 1 (e.g. only acidic aas):		return 1

	if ($ionisable_aas == 0) {
#		$pI = 999;
		$pI = 7;
	} elsif ($loop == 1) {
		$pI= 14;
	} elsif ($loop == 27) {
		$pI = 1;
	} else {
		$pI = ($pH+$pH_dec) + $pH_dec*($protein_charge/($protein_charge-$protein_charge_last));
	}

$pH		= 8.75;		# add a specific charge calculation for selected pH
$charge_here	= 0;
foreach $aa (@aa_keys) {
	if ($aanum{$aa} > 0) {
		$fac = $aafactor{$aa};
		$pH_pK = $pH - $aapK{$aa};
		if ($aa eq "X") {
			$charge_inc = -1.0 + $fac * (1/(1+10**($fac*$pH_pK)));
		} else {
			$charge_inc = $fac * (1/(1+10**($fac*$pH_pK)));
		}
		$charge_here += $aanum{$aa}*$charge_inc;
	}
}
# printf OUT "\n$pdbwords[0] $pH ";
# printf OUT "%7.3f\n", $charge_here;

	return $pI;
}

#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#----------------------
# SUBROUTINE protein_KD
#----------------------
# pass a list of single letter amino acid codes, return the pI (pH at
# neutral charge) for this protein

sub protein_KD {

#	KD hash is Kyte-Doolittle hydropathy index : JMB (1982) 157:105
#	for FOld Index, scale/norm following Uversky : Proteins (2000) 41:415, (-4.5,+4.5) to (0,1)
#	use KDnorm to get KDnorm_mean = sum(KDnorm) / winlen

$KD {'A'}	=  1.8;
$KD {'C'}	=  2.5;
$KD {'D'}	= -3.5;
$KD {'E'}	= -3.5;
$KD {'F'}	=  2.8;
$KD {'G'}	= -0.4;
$KD {'H'}	= -3.2;
$KD {'I'}	=  4.5;
$KD {'K'}	= -3.9;
$KD {'L'}	=  3.8;
$KD {'M'}	=  1.9;
$KD {'N'}	= -3.5;
$KD {'P'}	= -1.6;
$KD {'Q'}	= -3.5;
$KD {'R'}	= -4.5;
$KD {'S'}	= -0.8;
$KD {'T'}	= -0.7;
$KD {'V'}	=  4.2;
$KD {'W'}	= -0.9;
$KD {'Y'}	= -1.3;
$KD {'X'}	=  0.0;

$KD {'a'}	=  1.8;
$KD {'c'}	=  2.5;
$KD {'d'}	= -3.5;
$KD {'e'}	= -3.5;
$KD {'f'}	=  2.8;
$KD {'g'}	= -0.4;
$KD {'h'}	= -3.2;
$KD {'i'}	=  4.5;
$KD {'k'}	= -3.9;
$KD {'l'}	=  3.8;
$KD {'m'}	=  1.9;
$KD {'n'}	= -3.5;
$KD {'p'}	= -1.6;
$KD {'q'}	= -3.5;
$KD {'r'}	= -4.5;
$KD {'s'}	= -0.8;
$KD {'t'}	= -0.7;
$KD {'v'}	=  4.2;
$KD {'w'}	= -0.9;
$KD {'y'}	= -1.3;
$KD {'x'}	=  0.0;

# Normalise the KD hydropathy values.

  @KD_keys	= keys %KD;
  foreach $keyhere (@KD_keys) {
    $KDtmp		= $KD{$keyhere} + 4.5;
    $KDnorm{$keyhere}	= $KDtmp / 9.0;
#    printf BLAH "KDnorm for $keyhere is $KDnorm{$keyhere}\n"; 
  }

# put the passed sequence into the aa_seq array

  @aa_seq = @_;

  $KDnorm_tot	= 0;
  $n_KD_residues	= 0;
  foreach $aa (@aa_seq) {
    if (exists $KD{$aa}) {
      $KDnorm_tot	= $KDnorm_tot + $KDnorm{$aa};
      $n_KD_residues++;
    }
  }
  if ($n_KD_residues != 0) {
    $KDnorm_average	= $KDnorm_tot / $n_KD_residues;
  } else {
    $KDnorm_average	= 999;		# no residues, value unset, pick up in parent code
  }

return $KDnorm_average;

}

#-------------------------------------------------------------------------
#----------------------
# SUBROUTINE protein_Q
#----------------------
# pass a list of single letter amino acid codes, return the net charge per residue

# ?? his_charged ENV variable can be used to set his to +1 rather than zero, not used currently ??

sub protein_Q {

$chr {'A'}	=  0.0;
$chr {'C'}	=  0.0;
$chr {'D'}	= -1.0;
$chr {'E'}	= -1.0;
$chr {'F'}	=  0.0;
$chr {'G'}	=  0.0;
$chr {'H'}	=  0.0;
$chr {'I'}	=  0.0;
$chr {'K'}	=  1.0;
$chr {'L'}	=  0.0;
$chr {'M'}	=  0.0;
$chr {'N'}	=  0.0;
$chr {'P'}	=  0.0;
$chr {'Q'}	=  0.0;
$chr {'R'}	=  1.0;
$chr {'S'}	=  0.0;
$chr {'T'}	=  0.0;
$chr {'V'}	=  0.0;
$chr {'W'}	=  0.0;
$chr {'Y'}	=  0.0;
$chr {'X'}	=  0.0;

$chr {'a'}	=  0.0;
$chr {'c'}	=  0.0;
$chr {'d'}	= -1.0;
$chr {'e'}	= -1.0;
$chr {'f'}	=  0.0;
$chr {'g'}	=  0.0;
$chr {'h'}	=  0.0;
$chr {'i'}	=  0.0;
$chr {'k'}	=  1.0;
$chr {'l'}	=  0.0;
$chr {'m'}	=  0.0;
$chr {'n'}	=  0.0;
$chr {'p'}	=  0.0;
$chr {'q'}	=  0.0;
$chr {'r'}	=  1.0;
$chr {'s'}	=  0.0;
$chr {'t'}	=  0.0;
$chr {'v'}	=  0.0;
$chr {'w'}	=  0.0;
$chr {'y'}	=  0.0;
$chr {'x'}	=  0.0;

if ((exists $ENV{his_charged}) and ($ENV{his_charged} eq 'yes')) {
  $chr {'H'}	=  1.0;
  $chr {'h'}	=  1.0;
}

# put the passed sequence into the aa_seq array

  @aa_seq = @_;

  $Q_tot		= 0;
  $n_Q_residues		= 0;
  foreach $aa (@aa_seq) {
    if (exists $chr{$aa}) {
      $Q_tot	= $Q_tot + $chr{$aa};
      $n_Q_residues++;
    }
  }
  if ($n_Q_residues == 0) {
    $Q_mean		= 999;
  } else {
    $Q_mean		= $Q_tot / $n_Q_residues;
  }
  
  return $Q_mean;
}


#-------------------------------------------------------------------------
#----------------------------
# SUBROUTINE protein_disorder
#----------------------------
# pass a list of single letter amino acid codes, return the disorder calculation

sub protein_disorder {

# dis hash is Rune and Linding disorder propensity : GlobPlot : NAR (2003) 31:3701

$dis {'A'}	= -0.26154;
$dis {'C'}	= -0.01515;
$dis {'D'}	=  0.22763;
$dis {'E'}	= -0.20469;
$dis {'F'}	= -0.22557;
$dis {'G'}	=  0.43323;
$dis {'H'}	= -0.00122;
$dis {'I'}	= -0.42224;
$dis {'K'}	= -0.10009;
$dis {'L'}	= -0.33793;
$dis {'M'}	= -0.22590;
$dis {'N'}	=  0.22989;
$dis {'P'}	=  0.55232;
$dis {'Q'}	= -0.18768;
$dis {'R'}	= -0.17659;
$dis {'S'}	=  0.14288;
$dis {'T'}	=  0.00888;
$dis {'V'}	= -0.38618;
$dis {'W'}	= -0.24338;
$dis {'Y'}	= -0.20751;
$dis {'X'}	=  0.0;

$dis {'a'}	= -0.26154;
$dis {'c'}	= -0.01515;
$dis {'d'}	=  0.22763;
$dis {'e'}	= -0.20469;
$dis {'f'}	= -0.22557;
$dis {'g'}	=  0.43323;
$dis {'h'}	= -0.00122;
$dis {'i'}	= -0.42224;
$dis {'k'}	= -0.10009;
$dis {'l'}	= -0.33793;
$dis {'m'}	= -0.22590;
$dis {'n'}	=  0.22989;
$dis {'p'}	=  0.55232;
$dis {'q'}	= -0.18768;
$dis {'r'}	= -0.17659;
$dis {'s'}	=  0.14288;
$dis {'t'}	=  0.00888;
$dis {'v'}	= -0.38618;
$dis {'w'}	= -0.24338;
$dis {'y'}	= -0.20751;
$dis {'x'}	=  0.0;

# put the passed sequence into the aa_seq array

  @aa_seq = @_;

  $dis_tot	= 0;
  $ndis_tot	= 0;
  foreach $aa (@aa_seq) {
    if (exists $dis{$aa}) {
      $dis_tot += $dis{$aa};
      $ndis_tot++;
    }
  }

  if ($ndis_tot == 0) {
    $dis_mean		= 999;
  } else {
    $dis_mean		= $dis_tot / $ndis_tot;
  }
  
  return $dis_mean;
}


#-------------------------------------------------------------------------
#-----------------------
# SUBROUTINE protein_ent
#-----------------------
# pass a list of single letter amino acid codes, return the pI (pH at
# neutral charge) for this protein

sub protein_ent {

@aa_single_here	= ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");

# put the passed sequence into the aa_seq array

  @aa_seq = @_;

  for $aa (@aa_single_here) {
    $aa_count_here{$aa}	= 0;
  }

  $n_ent		= 0;
  foreach $aa (@aa_seq) {
    $AA_UC		= uc ($aa);			# should be upper case already, but only want to count 20 aas, not 40
    if (exists $aa_count_here{$AA_UC}) {
      $aa_count_here{$AA_UC}++;
      $n_ent++;
    }
  }

  if ($n_ent == 0) {
    $ent_sum		= 999;
  } else {
    $ent_sum		= 0;
    $frac_sum		= 0;
    for $aa (@aa_single_here) {
      if ($aa_count_here{$aa} != 0) {
        $frac_here	= $aa_count_here{$aa} / $n_ent;
        $incr_here	= -$frac_here * (log($frac_here) / log(2.0));
        $ent_sum		+= $incr_here;
        $frac_sum		+= $frac_here;
      }
    }
  }

  return $ent_sum;
}


#-------------------------------------------------------------------------
#------------------------------
# SUBROUTINE protein_betastrand
#------------------------------
# just returning the mean beta propensity - ?? could do more with 3 states

sub protein_betastrand {

# a very simple Chou-Fasman type ss prediction - see Costantini A et al (2006) BBRC 342:441-451

# read in the propensities

$sp92		= ' ';
for ($s=1; $s<=91; $s++) { $sp92 .= ' '; }
$sp81		= ' ';
for ($s=1; $s<=80; $s++) { $sp81 .= ' '; }

# no need for time_stamp - this is data, for read only
# use Cwd qw(getcwd);
# my $dcwd = getcwd;
use File::Basename;
my $dcwd = dirname(__FILE__);

$file_here		= $dcwd . '/ss_propensities.txt';
open (PROPENSITIES, "<$file_here") or die "cannot open $file_here\n";
$naa_here		= 0;
while ($line 		= <PROPENSITIES>) {
	chomp $line;
	my @words 	= split (" ",$line);
	if (exists $words[0]) {
	    if ($words[0] eq 'data') {			# amino acid data line
		$naa_here++;
		$aa		= $words[2];
		$p_alpha{$aa}	= $words[4];
		$p_beta{$aa}	= $words[5];
		$p_coil{$aa}	= $words[6];
	    }
	}
}
# printf BLAH "propensities read in for $naa_here amino acids\n";
close (PROPENSITIES);

# put the passed sequence into the aa_seq array

  @aa_seq = @_;

  $sum_alpha	= 0;
  $sum_beta	= 0;
  $sum_coil	= 0;
  $nseen	= 0;
  foreach $aa (@aa_seq) {
    if (exists $p_alpha{$aa}) {
      $sum_alpha	+= $p_alpha{$aa};
      $sum_beta		+= $p_beta{$aa};
      $sum_coil		+= $p_coil{$aa};
      $nseen++;
    }
  }
  if ($nseen == 0) {
    $beta_mean		= 999;
  } else {
    $beta_mean		= $sum_beta / $nseen;
  }

  return $beta_mean;
}

