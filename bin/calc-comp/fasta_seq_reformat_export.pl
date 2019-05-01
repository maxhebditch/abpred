#!/usr/bin/perl

$file_here		= 'blah.txt';
open (BLAH, ">>$file_here");		# file for comments, throughout server code

open (IN, "< reformat.in") or die "cannot open reformat_in\n";
@in=<IN>;
close (IN);

$file_here	= 'reformat_out';
open (OUT, ">$file_here") or die "cannot open reformat_out\n";

$seq_stored	= 'no';
$skip_seq	= 'no';
foreach $line (@in) {
  chomp $line;
  @words	= split (" ",$line);
  if ((substr ($line,0,1) eq '>') or ($words[0] eq '//')) {		# new sequence or a // terminator
    if ($seq_stored eq 'yes') {			# process and write out the previous sequence, if no errors
      $seqUC	= uc($seq_here);		# uppercase the sequence
      $seqlen	= length($seqUC);
      $seq_new	= "";
      $naa	= 0;
$n = 1;	# ?? VERY ODD does the loop not fix starting n at 1 ??
      for ($n==1; $n<=$seqlen; $n++) {
        $char_here	= substr($seqUC,$n-1,1);
        if (($char_here ne ' ') and ($char_here ne '*')) {		# ignore blanks and asterisks
	  if ($char_here =~ /[ACDEFGHIKLMNPQRSTVWY]/) {
	    $seq_new	.= $char_here;
	    $naa++;
#	  } else {		# not happy, but we can get <cr> at end, and easier to take only aa than skip allowed non-aa
#	    printf BLAH "**** stopping - fasta_seq_reformat - non 20 aa detected for protein ID = $current_id\n";
#	    $skip_seq	= 'yes';
	  }
	}
      }
      if ($skip_seq ne 'yes') {
        if ($naa >= 21) {
          printf OUT "$current_id\n";		# print the stored id
          printf OUT "$seq_new\n";
        } else {					# no aas
	  printf BLAH "**** stopping - fasta_seq_reformat - shorter than 21 aas for protein ID = $current_id\n";
        }
      }
      $skip_seq		= 'no';
    }
    if (substr ($line,0,1) eq '>') {
      $current_id	= $line;		# store the new seq ID - printer later, if seq is cleared
      $seq_stored	= 'yes';		# set up for new sequence
      $seq_here		= '';
    }
  } else {					# add on to curr seq, if not blank
    $len_here	= length($line);
    if ($len_here != 0) { 
      $seq_here .= $line;
    }
  }
}

if ($seq_stored eq 'yes') {			# print the last seq if it conforms
  $seqUC	= uc($seq_here);		# uppercase the sequence
  $seqlen	= length($seqUC);
  $seq_new	= "";
  $naa		= 0;
$n = 1;	# ?? VERY ODD does the loop not fix starting n at 1 ??
  for ($n==1; $n<=$seqlen; $n++) {
    $char_here	= substr($seqUC,$n-1,1);
    if (($char_here ne ' ') and ($char_here ne '*')) {		# ignore blanks and asterisks
      if ($char_here =~ /[ACDEFGHIKLMNPQRSTVWY]/) {
        $seq_new	.= $char_here;
	$naa++;
      }
    }
  }
  if ($skip_seq ne 'yes') {
    if ($naa >= 21) {
      printf OUT "$current_id\n";			# print the stored id
      printf OUT "$seq_new\n";
    } else {					# no aas
      printf BLAH "**** stopping - fasta_seq_reformat - shorter than 21 aas for protein ID = $current_id\n";
    }
  }
  $skip_seq	= 'no';
}

close (OUT);
close (BLAH);

exit;
