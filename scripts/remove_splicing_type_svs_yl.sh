# Detect and remove SVs that are actually splicing events
#
# Usage:
#   remove_splicing_type_svs.sh SVs.bedpe exons.bed
#
# Inputs:
#   SVs.bedpe       SVs in a BEDPE format
#   exons.bed       Exons in BED format. Fourth column is a unique ID for exons
#                   and the fifth column is a unique ID for genes. This file
#                   must have exactly six columns. 
# 
# Description:
#   Pseudogenes and cDNA both contain DNA with exon-exon junctions. These
#   junctions manifest as clusters of read pairs that support deletions
#   between exons of the affected genes. 
#
#   The script first identifies deletion-type SVs that are bridging
#   deletions. Genes with >= 2 exon-exon are assumed to be cDNA or pseudogene
#   insertion artefacts. Finally, all exon-exon bridging on cDNA/pseudogene
#   insertion artefact genes are removed. 
#
#   SV breakpoints are slopped for 10bp before overlapping with exons. 
#
#   'bedtools' must be installed and available.
#


BEDPE=$1
EXONS=$2
SAMPLE=${BEDPE##*/}
SAMPLE=${SAMPLE/.*}
PSEUDOGENE_INSERTION=$3

# Find overlaps between SV breakpoints and exons. Slop 50bps. 
awk '$1!~/^chrom/' $BEDPE \
|  perl -ane 'next if /^chrom/; print'   \
| perl -ane 'print join("\t", $F[0], $F[1]-50, $F[2]+50, $F[6], $F[7], $F[8], $F[3], $F[4]-50, $F[5]+50, $F[9]) . "\n"'  \
| bedtools intersect -a stdin -b $EXONS -loj  \
| perl -ane 'next if $F[10] eq "."; print'  \
| cut -f-16  \
| perl -ane 'print join("\t", @F[6..8, 3, 4, 9, 0..2, 5, 10..$#F]) . "\n"'  \
| bedtools intersect -a stdin -b $EXONS -loj  \
| perl -ane 'next if $F[16] eq "."; print'  \
| cut -f-22   \
| perl -ane 'print join("\t", @F[6..8, 0..2, 3, 4, 9, 5, 10..$#F]) . "\n"'  \
| perl -ane 'next unless $F[8] eq "+" and $F[9] eq "-"; print'  \
> $BEDPE.vs_exon_overlaps.txt

# Find genes with >= 2 exon-exon junctions. The SVs must be overlapping with distinct exons. 
detect_artefact_genes='
    while (<>) {
        chomp;
        @F = split /\t/;
        next if ($F[0] ne $F[3]  or  $F[8] ne "+"  or  $F[9] ne "-");  # Skip if SV is not deletion-type
        next if ($F[13] eq $F[19]);  # Skip if both ends overlap the same exon
        next if ($F[14] ne $F[20]);  # Skip if the exons are not on the same gene

        $exon_exon_counts_of_gene{$F[14]}++;
    }

    for (keys %exon_exon_counts_of_gene) {
        if ($exon_exon_counts_of_gene{$_} >= 2) {
            print "$_ $exon_exon_counts_of_gene{$_}\n";
        }
    }
'
perl -e "$detect_artefact_genes" $BEDPE.vs_exon_overlaps.txt  \
> $BEDPE.artefact_genes.txt

# Get pseudogene insertions of this sample
grep $SAMPLE $PSEUDOGENE_INSERTION \
| cut -f3  \
| perl -pe 'chomp; $_ = $_ . "\\.\n"'  \
> $BEDPE.inserted_pseudogenes.txt
grep -f $BEDPE.inserted_pseudogenes.txt $EXONS  \
| cut -f5  \
| sort | uniq  \
| perl -ne 'chomp; print "$_ 1000\n"'  \
> $BEDPE.artefact_genes.2.txt

cat $BEDPE.artefact_genes.txt $BEDPE.artefact_genes.2.txt \
| sort -u -k1,1  \
> $BEDPE.artefact_genes.3.txt

# Remove exon-exon deletion-type artefacts
remove_artefact_svs='
    $_ = `cut -f1 $ARGV[0]`;
    for (split /\s+/) {
        $is_bad_gene{$_} = 1;
    }
    $filter_all_genes = 0;
    if (scalar(keys %is_bad_gene) >= 3) {
        $filter_all_genes = 1;
    }

    open EXON_OVERLAPS, $ARGV[1] or die $!;
    %svs_to_be_skipped = ();
    while (<EXON_OVERLAPS>) {
        chomp;
        @F = split /\t/;
        
        # Updated 2016-09-11:
        # - From confirmed pseudogene insertions remove all SVs in exons. 
        if (exists($is_bad_gene{$F[14]}) or exists($is_bad_gene{$F[20]})) {
            $svs_to_be_skipped{$F[6]} = 1;  # Remove SV as long as it overlaps with an exon if the gene is confirmed pseudogene
        }
        elsif ($filter_all_genes) {
            next if ($F[0] ne $F[3]  or  $F[8] ne "+"  or  $F[9] ne "-");  # Skip if SV is not deletion-type
            next if ($F[13] eq $F[19]);  # Skip if both ends overlap the same exon
            next if ($F[14] ne $F[20]);  # Skip if the exons are not on the same gene
            $svs_to_be_skipped{$F[6]} = 1;
        }
    }
    close EXON_OVERLAPS;
    # print "$filter_all_genes\n";
    # print scalar(keys %is_bad_gene) . "\n";
    # print scalar(keys %svs_to_be_skipped) . "\n";

    open BEDPE, $ARGV[2] or die $!;
    while (<BEDPE>) {
        next if /^chrom/;
        @F = split /\t/;
        next if exists($svs_to_be_skipped{$F[6]});
        print;
    }
    close BEDPE;
'
perl -e "$remove_artefact_svs" $BEDPE.artefact_genes.3.txt $BEDPE.vs_exon_overlaps.txt $BEDPE  \
> ${BEDPE/.bedpe/.no_exon_exon_artefacts.bedpe}

cat $BEDPE.vs_exon_overlaps.txt | cut -f 7 | sort -u

