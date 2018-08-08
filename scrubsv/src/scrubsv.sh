#!/bin/bash

# $1 SV VCF
# $2 SNP VCF
phase-haplotypes () {
    svf="$1"
    svf_prefix=${svf%%.vcf.gz}
    svf_prefix=${svf_prefix%%.vcf}
    
    snpf="$2"

    hapcut_fn="${svf_prefix}.Hapcut.vcf"

    bam_f="$3"

    # output SNP VCF w/ phasing
    zcat <"$svf" |awk '!($10 ~ "0/3" || $10 ~ "1/4"|| $10 ~ "2/3"||$10 ~ "1/3")' > "$hapcut_fn"

    ##use extractHAIRS to convert BAM file to the compact fragment file format containing only haplotype-relevant information
    frag_fn="${svf_prefix}.fragment_file"

    extractHAIRS \
    --bam "${bam_f}" \
    --VCF "$hapcut_fn" --out "$frag_fn"

}

flag-overlap () {
    echo "This looks for overlaps between SVs"
}

segment-svs-vcf () {
    echo "This slices the svs vcf into manageable chunks for downstream analysis"
    # Output list of slice vcf files.
}

final-scrub () {
    echo "This takes as input "
    # May not need
}



main() {

    echo "Value of svs_vcf: '$svs_vcf'"
    echo "Value of snp_vcf: '$snp_vcf'"
    echo "Value of mapping_bam: '$mapping_bam'"

    dx download "$svs_vcf" -o "$svs_vcf_name"

    dx download "$snp_vcf" -o "$snp_vcf_name"

    dx download "$mapping_bam" -o "$mapping_bam_name"



    phased_vcf="$svs_vcf_name".phased.vcf
    phase-haplotypes "$svs_vcf_name" "$snp_vcf_name" > "$phased_vcf" &

    overlap_flags_tsv="$svs_vcf_name".overlap.flag.tsv
    flag-overlap


    flagged_vcf=$(dx upload flagged_vcf --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output flagged_vcf "$flagged_vcf" --class=file
}
