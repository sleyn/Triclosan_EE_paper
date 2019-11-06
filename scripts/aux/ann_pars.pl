#!/usr/bin/perl -w
use strict;

#Breaks variants in the same position in SNPEff vcf file to the separate strings

open ANNIN, "ann_1.vcf";
open ANNOUT, ">ann.vcf";

my @vcfs = glob("*AO*.vcf");
my %ann_vars;

foreach my $vcf_file ( @vcfs ){
    open VCF, $vcf_file;
    while(<VCF>){
		if( $_ =~ m/^#/ ){
		    next;
		}
		
		my @items = split /\t/, $_;
        my @alts = split /,/, $items[4];
        foreach my $alt ( @alts ){
            $ann_vars{"$items[0]\t$items[1]\t$items[3]\t$alt"} = 1;
        }
    }
}

my $ccc = 0;
my %check = ();    #check to exclude duplicates

while(<ANNIN>){
    $ccc++;
    chomp;
    if( $_ =~ m/^#/ ){
        print ANNOUT "$_\n";
    }else{
        my @items = split /\t/, $_;
        
        my @alts = split /,/, $items[4];            #alternative variants
        $items[7] =~ s/^(.+ANN=)//;
        my $prefix = $1;
        my @alt_anns = split /,/, $items[7];        #alternative annotations
        
        for( my $k = 0; $k < scalar(@alts); $k++){
            if( exists( $check{"$items[0]\t$items[1]\t$items[3]\t$alts[$k]"}) ){
                next;
            }
            
            if( exists( $ann_vars{"$items[0]\t$items[1]\t$items[3]\t$alts[$k]"} )){             #if variant was observed in the initial vcfs
                if( scalar(@alts) > 1 ){                                                        #if there is several alternatives  
                    my $string = $items[0];
                    for( my $j = 0; $j < scalar(@alt_anns); $j++ ){
                        if( exists( $check{"$items[0]\t$items[1]\t$items[3]\t$alts[$k]"}) ){
                            next;
                        }
                        if( $alt_anns[$j] =~ m/^$alts[$k]\|/ ){
                            for( my $i = 1; $i < scalar(@items); $i++){
                                if( $i == 4 ){
                                    $string .= "\t$alts[$k]";
                                    next;
                                }
                                if( $i == 7){
                                    
                                    $string .= "\t$prefix"."$alt_anns[$j]";
                                    next;
                                }
                                $string .= "\t$items[$i]";
                            }
                            print ANNOUT "$string\n";
                            $check{"$items[0]\t$items[1]\t$items[3]\t$alts[$k]"} = 1;
                        }
                    }
                }else{
                    print ANNOUT "$_\n";
                }
            }
        }
    }
}


