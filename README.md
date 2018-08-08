![alt text](ScrubSV_logo.png)

## ScrubSV
A toolkit to detect and flag potentially false SV's based on SNP's  and coverage.
Hackathon team: Lead: Fritz Sedlazeck - SysAdmin: Steve Osazuwa - Programmers: Priya Kritivasan, Kshithija Nagulapalli, David Oliver, Seungyeul Yoo

## How to use
Input requirements
1. bam file
2. Structural Variants in vcf format
3. Gold Standard SV callset 

## Installation

`git clone -r https://github.com/NCBI-Hackathons/ScrubSV.git`
`cd ScrubSV`

## Install HapCUT2
Refer to https://github.com/vibansal/HapCUT2 for installation instructions

## Install SURVIVOR
Refer to https://github.com/fritzsedlazeck/SURVIVOR for installation instructions. Move the executable to the current working directory.

## Dependencies: vcftools and Genomic Ranges
To install vcftools
`sudo apt-get vcftools`

To install GenomicRanges
`source("https://bioconductor.org/biocLite.R")`
`biocLite("GenomicRanges")`
`library(GenomicRanges)`

## Slides from the hackathon
We created and presented this method over the NYGC NCBI Hackathon in August 2018. 
[Here is the link to the slides](https://docs.google.com/presentation/d/16ZwtBfEyv7mvlIw1uxhgUtazK3lBlYRjU0Ys9OCQss8/edit?usp=sharing)



## Software Workflow Diagram
We included a workflow of our pipeline describing the modules in the repo. Feel free to check it out. 


