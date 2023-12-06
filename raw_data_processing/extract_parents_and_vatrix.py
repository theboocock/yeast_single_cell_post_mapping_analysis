# Copyright (c) 2021 Boocock James <james.boocock@otago.ac.nz>
# Author: Boocock James <james.boocock@otago.ac.nz>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import argparse

vatrix_cmd="/u/home/s/smilefre/project-kruglyak/sceqtl/2021/software/vartrix/target/release/vartrix"

def get_parent_list(parent_list):
    parent_dict = {} 
    with open(parent_list) as in_f:
        for line in in_f:
            l_s = line.strip().split(" ")
            p1 = l_s[0]
            p2 = l_s[1]
            parent = l_s[2]
            parent_dict[parent] = [p1,p2]
    return parent_dict

def make_sample_list(parent_dict, crosses):
    """
        
    """
    # Comma seperated list to list
    crosses = crosses.split(",")
    crosses = " ".join(crosses)
    parents = []
    if "all" in crosses:
        for key, value in parent_dict.items():
            parents.extend(value)
    else:
        for cross in crosses.split():
            parents.extend(parent_dict[cross])
    print(parents)
    parents = list(set(parents))
    return(parents)

import subprocess
__bcftools_cmd__="""bcftools view -v snps -m 2 -M 2 -s {samples} {vcf} | bcftools filter -i "MAC > 0" > parents.vcf"""

def run_bcftools_filter(vcf, parents):
    """
        Run bcftools filter
    """
    
    parents = ",".join(parents)
    bctools_cmd = __bcftools_cmd__.format(samples=parents,vcf=vcf)
    subprocess.check_call(bctools_cmd, shell=True)

import os
ref_fasta="/u/home/s/smilefre/project-kruglyak/sceqtl/2021/ref/yeast/fasta/genome.fa"
vatrix_cmd="""
/u/home/s/smilefre/project-kruglyak/sceqtl/2021/software/vartrix/target/release/vartrix --bam=possorted_genome_bam.bam \
    --cell-barcodes=filtered_feature_bc_matrix/barcodes.tsv \
    --out-matrix=alt_counts.mtx \
    --ref-matrix=ref_counts.mtx \
    --out-variants=out_var.txt \
    --fasta=/u/home/s/smilefre/project-kruglyak/sceqtl/2021/ref/yeast/fasta/genome.fa \
    --vcf=parents.vcf \
    --threads=16 \
    --scoring-method=coverage \
    --primary-alignments \
    --mapq=20 \
    --umi=TRUE \
    --log-level=info
"""

def run_vartrix():
    """
        Run VATRIX
    """
    if os.path.exists("filtered_feature_bc_matrix/barcodes.tsv.gz"):
        unzip_barcodes_commannd=""" gunzip filtered_feature_bc_matrix/barcodes.tsv.gz """
        subprocess.check_call(unzip_barcodes_commannd, shell=True)
    subprocess.check_call(vatrix_cmd,shell=True)

def main():
    parser=argparse.ArgumentParser(description="Extract parental genotypes from file")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--crosses",required=True)
    parser.add_argument("--parent-list",default="/u/home/s/smilefre/project-kruglyak/sceqtl/2021/ref/parent_list.txt")
    parser.add_argument("--threads",default=16)
    parser.add_argument("--cellranger-outdir")
    args=parser.parse_args()
    crosses = args.crosses
    vcf_input = os.path.abspath(args.vcf)
    parent_list = args.parent_list 
    parent_dict = get_parent_list(parent_list)
    parents = make_sample_list(parent_dict, crosses)
    cellranger_outdir = os.path.abspath(os.path.join(args.cellranger_outdir,"outs"))
    os.chdir(cellranger_outdir)
    run_bcftools_filter(vcf_input, parents)
    run_vartrix() 

if __name__=="__main__":
    main()
