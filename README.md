# SCRIPT FOR NGS (ChIP-seq) ANALYSIS
### Change file format from bam to .100.bw with parse2wig+

    build=hg38
    Dir=/work3/Database/Database_fromDocker/Referencedata_$build
    gt=$Ddir/genometable.txt
    gene=$Ddir/gtf_chrUCSC/chr.gene.refFlat
    sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/churros.0.13.2.sif"
    #definition of common variables
    gt=/work/Database/UCSC/mm39/genome_table
    post=.sort.bam  
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


### Peak calling with Drompa+


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
      #Definine directory of data

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

### Commmand to output Graphic Peak comparison between delta and ut sample using drompa+ (PFD file) exemple:

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
ls 5000 = changes bin size to 5kbp to show broad peaks


### Peak comparison with compare_bs, grouping remaining peaks after the depletion

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
                  
adding command "-and" output peaks of sample 1 that overlap to peaks of sample 2, so this is how peaks that remains after depletion can be grouped

### visualize results (exemple)
    $ parsecomparebs.pl peak/comparison_progmmousechipseq_delta_CTCF_CTCF.bed

this command outputs the statistics in a single line speced by tabs

### Grouping peaks: new and lost after deplation with compare_bs

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

compare_bs outputs the peaks that are in sample 1 but not in sample 2, if sample 1 is UT(untreated) and sample 2 is Delta(depletion) sample the output is group of lost peaks. To group new peaks is sufficent to exchange the order of the samples and create a new directory.           

### Peak Annotation 
Peak distribution using compare_bed2tss     #!/bin/bash

     build=mm39
     Ddir=/work3/Database/Database_fromDocker/Referencedata_$build
     gt=$Ddir/genometable.txt
     gene=$Ddir/gtf_chrUCSC/chr.gene.refFlat

     sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/churros.0.13.2.sif"

     #Definiion of marks and conditions
     conditions=("delta_CTCF" "delta_WAPL" "delta_RAD21")
     marks=("CTCF" "SMC1A" "SMC3" "RAD21")

     #Cicle for 
     for condition in "${conditions[@]}"; do
                  for mark in "${marks[@]}"; do

     $sing compare_bed2tss -g $gene --refFlat --gt $gt -m 2 -i \
          --bed peakcomparisonlost/comparison_progmmousechipseq_delta_CTCF_CTCF.bed \
          > peaklostdist/comparison_delta_CTCF_CTCF.bed

     done
     done

Creating heatmap with drompaplus

     sing="singularity exec --bind /work,/work2,/work3 /work3/SingularityImages/churros.0.13.2.sif"
     $sing drompa.heatmap.py \
    -o heatmap-aroundTSS \
     profile/aroundTES.PROFILE.averaged.ChIPread.ut_deltaCTCF_RAD21.tsv \
     profile/aroundTES.PROFILE.averaged.ChIPread.ut_deltaCTCF_SMC3.tsv \
     profile/aroundTSS.PROFILE.averaged.ChIPread.ut_deltaCTCF_RAD21.tsv \
     profile/aroundTSS.PROFILE.averaged.ChIPread.ut_deltaCTCF_SMC3.tsv

#-m 2 to put enerichement against genomic average, -i splits in exons and introns

### Rscript for pie charts of peak distribution 

      #Load necessary library
      library(ggplot2)
      library(gridExtra)

        #Data for delta under CTCF
        data_delta_ctcf <- list(
        CTCF = c(9772, 6336, 25988, 15768),
        SMC1A = c(5313, 3551, 11175, 9216),
        SMC3 = c(3176, 2321, 6671, 5435),
        RAD21 = c(6122, 4067, 12952, 9788)
        )

        #Data for delta under WAPL
        data_delta_wapl <- list(
        CTCF = c(2333, 1452, 6634, 4334),
        SMC1A = c(5454, 2448, 9198, 7857),
        SMC3 = c(1301, 648, 1995, 1101),
        RAD21 = c(909, 517, 1993, 1266)
        )

        #Data for delta under RAD21
        data_delta_rad21 <- list(
        CTCF = c(3527, 2350, 10446, 8112),
        SMC1A = c(6445, 4135, 13014, 11045),
        SMC3 = c(6556, 4594, 13946, 10697),
        RAD21 = c(5031, 3552, 12315, 9414)
          )

        #Labels for the categories
        labels <- c('Upstream', 'Downstream', 'Genic', 'Intergenic')
        colors <- c('#ff9999','#66b3ff','#99ff99','#ffcc99')

        #Function to create pie chart
         create_pie_chart <- function(data, protein_name, condition) {
            df <- data.frame(
          Category = labels,
          Values = data
         )
        ggplot(df, aes(x="", y=Values, fill=Category)) +
          geom_bar(width = 1, stat = "identity") +
          coord_polar("y", start=0) +
          scale_fill_manual(values = colors) +
          theme_void() +
          ggtitle(paste(condition, "-", protein_name)) +
          theme(plot.title = element_text(hjust = 0.5))
          }

      #Create pie charts for each condition
      plots_ctcf <- lapply(names(data_delta_ctcf), function(protein) create_pie_chart(data_delta_ctcf[[protein]], protein, "delta CTCF"))
      plots_wapl <- lapply(names(data_delta_wapl), function(protein) create_pie_chart(data_delta_wapl[[protein]], protein, "delta WAPL"))
      plots_rad21 <- lapply(names(data_delta_rad21), function(protein) create_pie_chart(data_delta_rad21[[protein]], protein, "delta RAD21"))

      # Arrange and save the plots
      grid.arrange(grobs = plots_ctcf, ncol = 2, nrow = 2)
      ggsave("delta_CTCF_pie_charts.png", arrangeGrob(grobs = plots_ctcf, ncol = 2, nrow = 2), width = 16, height = 8)

      grid.arrange(grobs = plots_wapl, ncol = 2, nrow = 2)
      ggsave("delta_WAPL_pie_charts.png", arrangeGrob(grobs = plots_wapl, ncol = 2, nrow = 2), width = 16, height = 8)

      grid.arrange(grobs = plots_rad21, ncol = 2, nrow = 2)
      ggsave("delta_RAD21_pie_charts.png", arrangeGrob(grobs = plots_rad21, ncol = 2, nrow = 2), width = 16, height = 8)


