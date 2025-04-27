#  This directory holds the base files for PoolingScreen
  
#  Future updates are likley to contain config files for easier use of the pipeline
  
  
#  The recommended computing resources would be to use a high resource computing cluster. 16 cores and 300 GB of ram are a god starting point as the pipeline includes assembly. 
#    Depending on the size of the raw sequencing dataset, this may influence the computing resources required. 
  
  
  
#  Main softwares and dependencies for the PoolingScreen pipeline and their versions used in the current version of the scripts:
  
        1. GNU awk: v4.0.2
        2. SRA toolkit: v2.10.9
        3. FastQC: v0.11.9
        4. MultiQC: v1.12
        5. Trimmomatic: v0.39
        6. SPAdes: v3.15.2
        7. DIAMOND: v2.0.15
        8. SeqKit: v0.10.0
        9. Krona: v2.8.1
        10. CheckV: v0.9.0
        11. TaxonKit: v0.3.0
  
  
#  Downloading of the NCBI databases for DIAMOND
  
#  Use the "virus_hermetia_databasecreationTAXID2022_05_15.sh" to download the full NCBI NR or NT databases and to also extract the sequences for the desired taxa using TaxonKit.
  
#  You can use the virus NCBI taxID "10239" and also the taxID of the host, if available and desired 
#  	Including the host will increase the time needed for the first DIAMOND BLASTx to run.
  
  
