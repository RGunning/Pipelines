#!/usr/bin/env perl

use strict;
my $ens_locs = {};

open (ENS_LOC, "<ens_genes_location.csv") or die "cant find location";
open (GENES,">ens_genes_and_promoter_distal.bed");
#open (PROMOTERS,">ens_promoters_only.csv");

while(<ENS_LOC>){
	chomp;
	# ENSMUSG00000064372,15356,15422,MT,-1,mt-Tp
	my ($ENS_ID,$START,$END,$CHR,$STRAND,$MARKER) = split /,/;
	$ens_locs->{$ENS_ID} = {
		CHR => $CHR,
		START => $START,
		END => $END,
		STRAND => $STRAND,
        MARKER => $MARKER
	};
}

print STDERR "recovered ".scalar(keys %$ens_locs)." ens ids\n";


foreach my $ens_id (keys %$ens_locs){
	my $ens_gene = $ens_locs->{$ens_id};
	my $marker = $ens_gene->{MARKER};
	my $chr= $ens_gene->{CHR};
	my $start= $ens_gene->{START};
	my $strand = $ens_gene->{STRAND};
	my $end= $ens_gene->{END};


# 	if($strand == "1"){
# 		$start = $start - 2000;
# 		if ($start < 0){$start = 0;}
# 	}else{
# 		$end = $end + 2000;
# 	}

    print GENES "chr${chr}\t${start}\t${end}\t${marker}\n";
#     print GENES "chr${chr}\t${promoter_start}\t${promoter_end}\t${marker}_Promoter\n";
#     print GENES "chr${chr}\t${distal_start}\t${distal_end}\t${marker}_Distal\n";
#    print GENES "chr${chr}\t${start}\t${end}\tGene\n";
#    print GENES "chr${chr}\t${promoter_start}\t${promoter_end}\tPromoter\n";
#    print GENES "chr${chr}\t${distal_start}\t${distal_end}\tDistal\n";
    #	print PROMOTERS "$ens_id,$marker,$chr,$promoter_start,$promoter_end,$strand\n";
}
