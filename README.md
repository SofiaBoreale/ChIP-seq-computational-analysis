# Peak calling with Drompa+


#!/bin/bash
#Creation of output directory
mkdir -p peak
DBdir=/work/Database/UCSC/mm39
gene=$DBdir/refFlat.NM.uniqTSS.txt

sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/churros.0.13.2.sif"


#Definition of common variables
gt=/work/Database/UCSC/mm39/genome_table
post=.100.bw  # $dir のファイルのpostfix
dir=/work2/CohesinProject/GSE178982_Hsieh_NatGenet2022/ChIP-seq/parse2wigdir+
#Definisciil percorso della directory dei dati

#Definitions of marks and conditions
conditions=("CTCF" "SMC1A" "SMC3" "RAD21")
marks=("delta_CTCF_IAA" "delta_CTCF_UT" "delta_WAPL_IAA" "delta_WAPL_UT" "delta_RAD21_IAA" "delta_RAD21_UT")

#Cicle for to create the variables and execute the command with drompa+
for condition in "${conditions[@]}"; do
          for mark in "${marks[@]}"; do
                      s="-i $dir/${mark}_${condition}ChIP$post,$dir/${mark}_CTCFinput$post,${mark}_${condition}"
                          echo "Executing drompa+ per ${condition} in condition ${mark}..."

                              #Execute the command
                                  $sing drompa+ PC_SHARP -g $gene $s -o peak/progmouse_chipseq_${condition}_${mark} --gt $gt --callpeak
                         done
                         done
