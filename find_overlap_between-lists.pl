#!/usr/bin/perl -w
#Script written by Wouter van den Berg
#Script finds C. elegans orthologs to a give list of C. briggsae genes, using the file - briggsae-elegans-orthologs-minimum-list-Sept-2018.txt

#Last modified on June 10 2023
#
#Comments
#-------------


#Declare files to be used in the script
$tagIN1 				= "input.txt";
$tagIN2 				= "briggsae-elegans-orthologs-list.txt";
$tagOUT 				= "output.txt";



open TAGIN1, $tagIN1 or die "Couldn't open the file containing tags: $!\n"; 
open TAGOUT, ">>$tagOUT" or die "Couldn't open the file tags: $!\n";
print TAGOUT "Cbr-WBGene\tCel-WBgene\tCel-common\Cbr_log2\n";

while (<TAGIN1>) {
	chomp $_;
	
	if(/#/) {print "Skipping the comment line in file 1.\n";} 
	else {
		@element1 	= split('\t', $_);
		$cbr_wb1 = $element1[0];
		$cbr_wb1   =~s/\s//g; 
		 $cbr_log2 = $element1[3];
		 $cbr_log2   =~s/\s//g; 


		open TAGIN2, $tagIN2 or die "Couldn't open the file containing tags: $!\n"; 



		while (<TAGIN2>) {
			chomp $_;
			
			if(/#/) {print "Skipping the comment line in file 2.\n";} 
			else {
				@element2 	= split('\t', $_); 
				$cbr_wb2 = $element2[1];
				$cbr_wb2 =~s/\s//g;
				$cel_wb2 = $element2[3];
				$cel_wb2 =~s/\s//g; 
				$cel_common2 = $element2[2];
				$cel_common2 =~s/\s//g; 
		
				
				if ($cbr_wb1 eq $cbr_wb2) {

						print "$cbr_wb1\t$cel_wb2\t$cel_common2\t$cbr_log2\n";
						print TAGOUT "$cbr_wb1\t$cel_wb2\t$cel_common2\t$cbr_log2\n";
				}
			}
		}			
	 }	
}
				
	
		