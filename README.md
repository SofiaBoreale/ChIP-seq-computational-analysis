# Change file format from bam to .100.bw with parse2wig+

build=hg38

Ddir=/work3/Database/Database_fromDocker/Referencedata_$build

gt=$Ddir/genometable.txt

gene=$Ddir/gtf_chrUCSC/chr.gene.refFlat

sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/churros.0.13.2.sif"

#definition of common variables
    
    gt=/work/Database/UCSC/mm39/genome_table
    post=.sort.bam  # $dir のファイルのpostfix
    dir=/work2/CohesinProject/GSE178982_Hsieh_NatGenet2022/ChIP-seq/bam #Definisci il percorso della directory dei dati
    

#definition of marks and conditions
      
      conditions=("CTCF" "SMC1A" "SMC3" "RAD21")
      marks=("delta_CTCF_IAA" "delta_CTCF_UT" "delta_WAPL_IAA" "delta_WAPL_UT" "delta_RAD21_IAA" "delta_RAD21_UT")

#cicle for to define variables and execute with drompa+

     for condition in "${conditions[@]}"; do
         for mark in "${marks[@]}"; do
                          $sing parse2wig+ -i $dir/${mark}_${condition}_input$post  -o ${mark}_${condition}input --gt $gt -p 8
                    done
                    done

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
                          

            $sing drompa+ PC_SHARP -g $gene $s -o peak/progmouse_chipseq_${condition}_${mark} --gt $gt --callpeak
                         done
                         done

# Commmand to output Graphic Peak comparison between delta and ut sample using drompa+ (PFD file) exemple:


#!/bin/bash

 mkdir -p drompa3

DBdir=/work/Database/UCSC/mm39

gene=$DBdir/refFlat.NM.uniqTSS.txt

sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/churros.0.13.2.sif"

gt=/work/Database/UCSC/mm39/genome_table


dir=/work2/CohesinProject/GSE178982_Hsieh_NatGenet2022/ChIP-seq/parse2wigdir+

        $sing drompa+ PC_SHARP -g $gene -i $dir/delta_CTCF_IAA_CTCFChIP.100.bw,$dir/delta_CTCF_IAA_CTCFinput.100.bw,deltaCTCF_CTCF \
         -i $dir/delta_CTCF_UT_CTCFChIP.100.bw,$dir/delta_CTCF_UT_CTCFinput.100.bw,utCTCF_CTCF -o drompa3 --gt $gt 
       --showitag 2 
       --ls 5000   
#ls 5000 = changes bin size to 5kbp to show broad peaks

#shows the input lines reads

#outputs in drompa3 directory

# Peak comparison with compare_bs, grouping remaining peaks after the depletion

#!/bin/bash

sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/churros.0.13.2.sif"

mkdir -p peakcomparison

 conditions=("delta_CTCF" "delta_WAPL" "delta_RAD21")
 marks=("CTCF" "SMC1A" "SMC3" "RAD21")

#Cicle for 

       for condition in "${conditions[@]}"; do
                  for mark in "${marks[@]}"; do
                                    
        $sing compare_bs -and -1 peak/progmouse_chipseq_${mark}_${condition}_UT.${condition}_UT_${mark}.peak.bed \
                 -2 peak/progmouse_chipseq_${mark}_${condition}_IAA.${condition}_IAA_${mark}.peak.bed \
                  > peakcomparison/comparison_progmmousechipseq_${condition}_${mark}.bed
                  done         
                  done
                  
#adding command -and i needed to output peaks of sample 1 that overlap to peaks of sample 2, so this is how peaks that remains after depletion can be grouped

# visualize results (exemple)
    $ parsecomparebs.pl peak/comparison_progmmousechipseq_delta_CTCF_CTCF.bed

#this command outputs the statistics in a single line speced by tabs

# Grouping peaks: new and lost after deplation with compare_bs

#!/bin/bash

sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/churros.0.13.2.sif"

#creates the output drectory

mkdir -p peakcomparisonlost

conditions=("delta_CTCF" "delta_WAPL" "delta_RAD21")

marks=("CTCF" "SMC1A" "SMC3" "RAD21")

    # Cicle for 
     for condition in "${conditions[@]}"; do
                  for mark in "${marks[@]}"; do
                                
        $sing compare_bs -1 peak/progmouse_chipseq_${mark}_${condition}_UT.${condition}_UT_${mark}.peak.bed \
                 -2 peak/progmouse_chipseq_${mark}_${condition}_IAA.${condition}_IAA_${mark}.peak.bed \
                  > peakcomparisonlost/comparison_progmmousechipseq_${condition}_${mark}.bed
                  done         
                  done

#compare_bs outputs the peaks that are in sample 1 but not in sample 2, if sample 1 is UT(untreated) and sample 2 is Delta(depletion) sample the output is group of lost peaks. To group new peaks is sufficent to exchange the order of the samples and create a new directory           

