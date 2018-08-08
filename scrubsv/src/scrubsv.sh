#!/bin/bash


assemble-haplotypes () {
    TODO
    # output SNP VCF w/ phasing
}

check-overlap () {
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




    flagged_vcf=$(dx upload flagged_vcf --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output flagged_vcf "$flagged_vcf" --class=file
}
