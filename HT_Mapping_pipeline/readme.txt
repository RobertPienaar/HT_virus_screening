#  Here are the scripts for the high-throughput virus screening by mapping.
  
  
#  The recommended computing resources would be to use a high-resource computing cluster. 
#    Typically this pipeline is run using 24 cores and 100 GB of RAM. 
#    However, the size of the dataset and reference sequence will influence this. 
#    For instance, the visualisations can easily be run on a laptop with 16GB of RAM and 8 cores/threads. 
#    The final visualisation of the BAM files in Rstudio on a laptop may take some time (more than 30 minutes) to complete per reference sequence if trying to visualise more than 50 datasets together, depending on the size of the BAM files. 
  
  
#  For the mapping, the versions of the software used:
  
    1. fastp: v0.23.2
    2. Bowtie2: v2.4.2
    3. Samtools: v1.9
  
  
#  Rpackages and their versions used in the current version of the scripts:
  
#  Genome and basic annotation visualisation
  
      1. Rbase: v4.2.2 
  
#  Mapping visualisation
  
      1. Rsamtools: v2.14.0
      2. ComplexHeatmap: v2.14.0
      3. dplyr: v1.1.4
      4. GenomicRanges: v1.50.2
      5. GenomicAlignments: v1.34.1
      6. devtools: v2.4.5
  
