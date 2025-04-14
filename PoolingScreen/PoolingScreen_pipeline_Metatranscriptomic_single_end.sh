#!/bin/bash

eval "$(conda shell.bash hook)"

D_And_T="date --universal +%A__%d-%B-%Y__%H:%M:%S"

# Tip: Place all the raw reads in the ProjectPATH and in their own directory which is given a "Projectname". You can use this Project name in the variable "Bioprojectname" for the loop to know what folder to use.

#-------------------------------- Software to run ----------------------------------#
# 1
Run_SRAtoolkit_fasterq_dump=YES
# 2
Run_fastQC_and_MultiQC=YES
# 3
Run_trimmomatic=YES
# 4
assembler_used=spades # don't change this unless you want to used a different assembler, if this is the case, then also change "Run_rnaSPAdes" to "NO"
contig_name_prefix=NODE_ # don't change this unless you want to used a different assembler, if this is the case, then also change "Run_rnaSPAdes" to "NO"
Run_rnaSPAdes=YES
# 5
Change_contig_names=YES
# 6
Run_diamond_spades_virus=YES
# 7
Update_Krona_taxonomy=NO
Run_Krona_on_assembled_contigs=YES
Run_Krona_project_on_assembled_contigs=YES
# 8
Extract_assembled_contigs_with_viral_hits=YES
# 9
Run_diamond_virus_hits_full_db=YES
# 10
Extract_viral_candidate_contigs=YES
# 11
Run_Krona_on_viral_candidates=YES
Run_Krona_project_on_viral_candidate_contigs=YES
# 12
Do_you_want_to_concatenate_fasta_files_for_CheckV=YES
# 13
Run_CheckV=YES
# 14
Run_Krona_project_on_viral_candidate_contigs_of_all_projects=YES

#-------------------------- Location of Working directory --------------------------#
WorkingPATH=/home/Dataanalysis/Curations/Metatranscriptomic_SE/Date
ProjectPATH=${WorkingPATH}/raw_reads
BioprojectName=(Project1 Project2 Project3) # You can place multiple bioproject names in the brackets, but just separate them by a single space. e.g. Project1 Project2

Date_of_pipeline_run="date +%d_%B_%Y"

#-------------------------- Main directories to make -------------------------------#
SRA_files_output_directory=/home/Dataanalysis/raw_transcriptome_reads
mkdir ${SRA_files_output_directory}

Pipeline_output_directory=${WorkingPATH}/Virus_Pipeline_`${Date_of_pipeline_run}`
mkdir ${Pipeline_output_directory} 

FastQC_output_directory=${Pipeline_output_directory}/FASTQC
mkdir ${FastQC_output_directory}

MultiQC_output_directory=${Pipeline_output_directory}/MultiQC
mkdir ${MultiQC_output_directory}

Trimmomatic_output_directory=${Pipeline_output_directory}/trimmomatic # Do not change this unless you know that your version accepts a directory name other than "trimmomatic"
mkdir ${Trimmomatic_output_directory}

rnaSPAdes_output_directory=${Pipeline_output_directory}/rnaSPAdes
mkdir ${rnaSPAdes_output_directory}

Extracted_output_contigs=${Pipeline_output_directory}/Extracted_outputs
mkdir ${Extracted_output_contigs}

renamed_assembly_outputs=${Extracted_output_contigs}/${assembler_used}_nodes
mkdir ${renamed_assembly_outputs}

Diamond_output_directory=${Extracted_output_contigs}/spades_diamond_outputs
mkdir ${Diamond_output_directory}

Krona_output_directory=${Extracted_output_contigs}/spades_diamond_kronaout
mkdir ${Krona_output_directory}

CheckV_output_directory=${Extracted_output_contigs}/spades_CheckV_outputs
mkdir ${CheckV_output_directory}

#---------------------------------- Pipeline Logfiles ------------------------------#

Pipeline_Logfile=${Pipeline_output_directory}/Virus_pipeline_`${Date_of_pipeline_run}`.logfile
Commands_check_log=${Pipeline_output_directory}/Virus_pipeline_commands_`${Date_of_pipeline_run}`.logfile

#------------------------------ Assembler parameters -------------------------------#

#When changing the file names of the output contigs of the assembler. e.g. rnaSPAdes will label the file "transcripts.fasta"
assembler_fasta_PREFIX=transcripts

InputPath_rnaSPAdes=${Trimmomatic_output_directory}  # Change this only if you do not want to use Trimmomatic.

# You can change the extensions for the inputs of rnaSPAdes here to your preference. If you are using Trimmomatic, the pipeline will automatically adjust the extensions of the Trimmomatic outputs to fit your choice.

InputFileExtension_rnaSPAdes=.fastq.qtrim.fastq #Forward reads or reaction 1 reads - change "fastq" to "fasta" if you are using "fasta" inputs for Trimmomatic
   
#___________________________________________________________________________________#

#----------------------------- Trimmomatic parameters ------------------------------#

Path_to_adapters=/home/Dataanalysis/Trimmomatic_Adapters
Adapters_file=${Path_to_adapters}/TruSeq2_3-SE.fa 

# list of adapter files
#TruSeq2_3 has the adapters from TruSeq2 and TruSeq3

##                  Nextera                   ##
# NexteraPE-PE.fa                              #
#::::::::::::::::::::::::::::::::::::::::::::::#

##                   TruSeq                   ##
#:::::::::::::::::#          #:::::::::::::::::#
#   Single-end    #          #   Paired-end    #
#:::::::::::::::::#          #:::::::::::::::::#
# TruSeq2-SE.fa   #          # TruSeq2-PE.fa   #
# TruSeq3-SE.fa   #          # TruSeq3-PE.fa   #
# TruSeq2_3-SE.fa #          # TruSeq3-PE-2.fa #
###################          # TruSeq2_3-PE.fa #
                             ###################
#::::::::::::::::::::::::::::::::::::::::::::::#

# You can edit this in the pipeline to your preference but this should be fine for RNAseq data - Based off Trinity v2.11.0
# ILLUMINACLIP:${Adapters_file}:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 #Don't remove the "#" at the front


#------------------------------- Diamond parameters --------------------------------#

# where the database files are stored, if you are using the full nr.dmnd database use 
Virus_DatabaseForDIAMOND=/home/Dataanalysis/databases/diamond/nr_virus_hermetia_db_2022_05_15.dmnd
Full_DatabaseForDIAMOND=/home/Dataanalysis/databases/diamond/nr_2022_05_15.dmnd
virus_database_shortname=hermetia_virusrndb
Full_database_shortname=Full_NR_db
Project_collection=BSF_SE_transcriptomes_by_date

# The extension for the input sample files, are they .fastq/.fq or .fasta/.fa? keep the period punctuation mark "." in the variable.
#--------------InputFileExtension=.fasta
Diamond_InputFileExtension=.fasta

# Approximate amount of memory which you would want DIAMOND to use. Remember, try to leave some memory available that you had originally requested as this value below is just an estimate, DIAMOND may use more or less.
DiamondMEMORY=300

### Do not change this unless it is causing an error, if it is causing an error, or to check, put a # in front of the "DiamondBlocksize=..." and just go to the diamond command and replace "${DiamondBlocksize}" and replace it with the number of GB that you wish to use for DIAMOND (remember to take 10% of your requested amount of RAM/memory and devide that number by 6, as Diamond will use roughly 6 times the amount that is placed in the "--block-size <number>" parameter.
DiamondBlocksize=$(echo "$((${DiamondMEMORY} / 6))")

#------------------------ Last few important parameters to change ---------------------------#

Taxa_accessions_list=/home/Dataanalysis/databases/taxonkit/downloaded/generated_dbs/nr_virus_id.acc.txt.txt # This is the file with the viral proteins accessions that you want to extract from the assemblies 

# Conda environments
Krona_environment=/home/environment-conda/kraken2
CheckV_environment=/home/environment-conda/checkv
CheckV_DB=/home/Dataanalysis/databases/checkv/checkv-db-v1.2

# Number ot threads/CPUs to be used
Num_threads=92
Amount_of_RAM=500
Format_type_FASTQC=fastq
InputFileExtension_FASTQC_and_Trimmomatic=fastq

#------------------ Software variables (how to call each program) ------------------# 
FasterqdumpPROG=fasterq-dump #version 2.10.9
FastQCPROG=fastqc #version 0.11.9
MultiQCPROG=multiqc #version 1.12
trimmomaticPROG=trimmomatic #version 0.39
SPAdesPROG=/home/Software/SPAdes_3.15.2_Linux/bin/spades.py #version 3.15.2
diamondProg=/home/Software/diamond/diamond #version 2.0.15
seqkitProg=seqkit #version 0.10.0
KronaBLASTImportPROG=ktImportBLAST #version 2.8.1
Krona_taxa_PROG=ktUpdateTaxonomy.sh #version 2.8.1
CheckV_PROG=checkv #version 0.9.0
#awk # GNU awk version 4.0.2

#============================================================================ Pipeline start ============================================================================#

echo "#____________ Virus Discovery Pipeline Log ____________#" > ${Pipeline_Logfile}
${D_And_T} >> ${Pipeline_Logfile}
echo "This logfile is to help you track the steps of the pipeline" >> ${Pipeline_Logfile}
echo "#____________ Virus Discovery Pipeline Command check Log ____________#" > ${Commands_check_log}
${D_And_T} >> ${Commands_check_log}
echo "This logfile is to help you troubleshoot any commands that may have been entered incorrectly" >> ${Commands_check_log}


#____________________________ Software versions ____________________________#
echo "#_______________________________ Software versions ______________________________#" >> ${Pipeline_Logfile}


if [ ${Run_SRAtoolkit_fasterq_dump} = YES ]
then
   echo "# SRAtoolkit fasterq dump" >> ${Pipeline_Logfile}
   ${FasterqdumpPROG} --version >> ${Pipeline_Logfile}
else
   echo "Fasterq dump not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Run_fastQC_and_MultiQC} = YES ]
then
   echo "# QC" >> ${Pipeline_Logfile}
   ${FastQCPROG} --version >> ${Pipeline_Logfile}
   ${MultiQCPROG} --version >> ${Pipeline_Logfile}
else
   echo "FastQC qnd MultiQC not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Run_trimmomatic} = YES ]
then
   echo "# Trimmomatic" >> ${Pipeline_Logfile}
   ${trimmomaticPROG} -version >> ${Pipeline_Logfile}
else
   echo "Trimmomatic not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Run_rnaSPAdes} = YES ]
then 
   echo "# Assembly in this run" >> ${Pipeline_Logfile}
   ${SPAdesPROG} --version >> ${Pipeline_Logfile}
else
   echo "Assembly not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Change_contig_names} = YES ]
then 
   echo "# contig names wil be changed to include project name and sample name in this run" >> ${Pipeline_Logfile}
else
   echo "contig names wil be changed to include project name and sample name not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Run_diamond_spades_virus} = YES ]
then 
   echo "# Diamond BLASTx using virus database on SPAdes output in this run" >> ${Pipeline_Logfile}
   ${diamondProg} --version >> ${Pipeline_Logfile}
else
   echo "Diamond BLASTx on SPAdes output using virus database not included in this run" >> ${Pipeline_Logfile}
fi 

if [ ${Extract_assembled_contigs_with_viral_hits} = YES ]
then 
   echo "# Assembled contigs with virus hits will be extracted in this run" >> ${Pipeline_Logfile}
   ${seqkitProg} version >> ${Pipeline_Logfile}
   awk --version >> ${Pipeline_Logfile}
else
   echo "Assembled contigs with virus hits extraction not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Run_Krona_on_assembled_contigs} = YES ]
then 
   echo "# Krona on assembled contigs in this run" >> ${Pipeline_Logfile}
   echo "check version manually" >> ${Pipeline_Logfile}
else
   echo "Krona on assembled contigs not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Run_Krona_project_on_assembled_contigs} = YES ]
then 
   echo "# Krona multi-sample on assembled contigs in this run" >> ${Pipeline_Logfile}
   echo "check version manually" >> ${Pipeline_Logfile}
else
   echo "Krona multi-sample on assembled contigs not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Run_diamond_virus_hits_full_db} = YES ]
then 
   echo "# Diamond BLASTx using full database on extracted viral hits output in this run" >> ${Pipeline_Logfile}
   ${diamondProg} --version >> ${Pipeline_Logfile}
else
   echo "Diamond BLASTx using full database on extracted viral hits output not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Extract_viral_candidate_contigs} = YES ]
then 
   echo "# Viral candidate contigs will be extracted in this run" >> ${Pipeline_Logfile}
   ${seqkitProg} version >> ${Pipeline_Logfile}
   awk --version >> ${Pipeline_Logfile}
else
   echo "Viral candidate contigs extraction not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Run_Krona_on_viral_candidates} = YES ]
then 
   echo "# Krona on viral candidates and full db contigs" >> ${Pipeline_Logfile}
   echo "check version manually" >> ${Pipeline_Logfile}
else
   echo "Krona on viral candidates and full db contigs not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Run_Krona_project_on_viral_candidate_contigs} = YES ]
then 
   echo "# Krona multi-sample on viral candidates and full db contigs" >> ${Pipeline_Logfile}
   echo "check version manually" >> ${Pipeline_Logfile}
else
   echo "Krona multi-sample on viral candidates and full db contigs not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Run_CheckV} = YES ]
then 
   echo "# CheckV on virus candidate contigs in this run" >> ${Pipeline_Logfile}
   echo "check version manually" >> ${Pipeline_Logfile}
else
   echo "CheckV not included in this run" >> ${Pipeline_Logfile}
fi

if [ ${Run_Krona_project_on_viral_candidate_contigs_of_all_projects} = YES ]
then 
   echo "# Krona multi-sample on viral candidates for all projects" >> ${Pipeline_Logfile}
   echo "check version manually" >> ${Pipeline_Logfile}
else
   echo "Krona multi-sample on viral candidates for all projects not included in this run" >> ${Pipeline_Logfile}
fi

echo "#________________________________________________________________________________#" >> ${Pipeline_Logfile}
#___________________________________________________________________________#

#Updating Krona Taxonomy

if [ ${Update_Krona_taxonomy} = YES ]
then
   echo "source activate ${Krona_environment}" >> ${Commands_check_log}
   source activate ${Krona_environment}
   echo "environment activated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   ${Krona_taxa_PROG} && >> ${Commands_check_log}
   echo "conda deactivate" >> ${Commands_check_log}
   conda deactivate
   echo "environment deactivated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
else
   echo "Krona taxonomy will not be updated in this run" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
fi

for ProjectName in ${BioprojectName[*]} 
do 
   InputPath1=${ProjectPATH}/${ProjectName}

#----------------------------------------------------------------------#
#            1. Converting SRA reads into fastq files                  #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#            1. Converting SRA reads into fastq files                  #\n#----------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   # fasterq dump
   if [ ${Run_SRAtoolkit_fasterq_dump} = YES ]
   then 
      echo "Script will run fasterq dump" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}

      mkdir ${SRA_files_output_directory}/${ProjectName}

      # To create a string (horizontal list) which has all the sample names with spaces inbetween each name 
      # if you want to do all of this in one line and not have any file remaining afterwards
      InputShortName=($(echo $(find ${InputPath1} -type f -name "*" -printf "%f ")))
      echo "input samples" >> ${Pipeline_Logfile}
      echo "${InputShortName}" >> ${Pipeline_Logfile}
      for Input in ${InputShortName[*]} 
      do 
         echo "###Input sample for fasterq dump is \"${Input}\"" >> ${Pipeline_Logfile}
         ##fasterq dump
         echo "##fasterq dump" >> ${Pipeline_Logfile}
         #fasterq dump command
         
         echo "${FasterqdumpPROG} ${InputPath1}/${Input} --outdir ${InputPath1}" >> ${Commands_check_log}
         ${FasterqdumpPROG} ${InputPath1}/${Input} --outdir ${InputPath1}
         mv ${InputPath1}/${Input} ${SRA_files_output_directory}/${ProjectName}/${Input}
         echo "fasterq dump complete for ${Input}" >> ${Pipeline_Logfile}
      done
   else
      echo "Script will NOT run fasterq dump" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi

#----------------------------------------------------------------------#
#            2. Assessing the quality of the sequencing                #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#            2. Assessing the quality of the sequencing                #\n#----------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   # FastQC
   if [ ${Run_fastQC_and_MultiQC} = YES ]
   then 
      echo "Script will run FastQC" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}

      mkdir ${FastQC_output_directory}/${ProjectName}

      # To create a string (horizontal list) which has all the sample names with spaces inbetween each name 
      # if you want to do all of this in one line and not have any file remaining afterwards
      InputShortName=($(sed -e "s+.${InputFileExtension_FASTQC_and_Trimmomatic}++g" <(echo $(find ${InputPath1} -type f -name "*.${InputFileExtension_FASTQC_and_Trimmomatic}" -printf "%f "))))
      echo "input samples" >> ${Pipeline_Logfile}
      echo "${InputShortName}" >> ${Pipeline_Logfile}
      for Input in ${InputShortName[*]} 
      do 
         echo "###Input sample for FastQC is \"${Input}\"" >> ${Pipeline_Logfile}
         ##FastQC
         echo "##FastQC" >> ${Pipeline_Logfile}
         #FastQC command
         
         echo "${FastQCPROG} --outdir ${FastQC_output_directory}/${ProjectName} --noextract --format ${Format_type_FASTQC} ${InputPath1}/${Input}.${InputFileExtension_FASTQC_and_Trimmomatic}" >> ${Commands_check_log}
         
         ${FastQCPROG} --outdir ${FastQC_output_directory}/${ProjectName} --noextract --format ${Format_type_FASTQC} ${InputPath1}/${Input}.${InputFileExtension_FASTQC_and_Trimmomatic}
         echo "FastQC complete for ${Input}" >> ${Pipeline_Logfile}
      done
      
      echo "Running MultiQC for FastQC outputs" >> ${Pipeline_Logfile}
      mkdir ${MultiQC_output_directory}/${ProjectName}
      
      echo "${MultiQCPROG} --outdir ${MultiQC_output_directory}/${ProjectName} --force --zip-data-dir --title ${ProjectName} ${FastQC_output_directory}/${ProjectName}" >> ${Commands_check_log}

      ${MultiQCPROG} --outdir ${MultiQC_output_directory}/${ProjectName} --force --zip-data-dir --title ${ProjectName} ${FastQC_output_directory}/${ProjectName}

   else
      echo "Script will NOT run FastQC and MultiQC" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi

#----------------------------------------------------------------------#
#                   3. QC and trimming of reads                        #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#                   3. QC and trimming of reads                        #\n#----------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   # Trimmomatic
   if [ ${Run_trimmomatic} = YES ]
   then 
      echo "Script will run Trimmomatic" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}

      mkdir ${Trimmomatic_output_directory}/${ProjectName}

      # To create a string (horizontal list) which has all the sample names with spaces inbetween each name 
      # if you want to do all of this in one line and not have any file remaining afterwards.
      InputShortName=($(sed -e "s+.${InputFileExtension_FASTQC_and_Trimmomatic}++g" <(echo $(find ${InputPath1} -type f -name "*.${InputFileExtension_FASTQC_and_Trimmomatic}" -printf "%f "))))
      echo "input samples" >> ${Pipeline_Logfile}
      echo "${InputShortName}" >> ${Pipeline_Logfile}
      for Input in ${InputShortName[*]} 
      do 
         echo "###Input sample for Trimmomatic is \"${Input}\"" >> ${Pipeline_Logfile}
         ${D_And_T} >> ${Pipeline_Logfile}
         ##Trimmomatic
         echo "##Trimmomatic" >> ${Pipeline_Logfile}
         #Trimmomatic command

         echo "${trimmomaticPROG} SE -threads ${Num_threads} -trimlog ${Trimmomatic_output_directory}/${Input}.logfile -summary ${Trimmomatic_output_directory}/${Input}_statssummary.txt ${InputPath1}/${Input}.${InputFileExtension_FASTQC_and_Trimmomatic} ${Trimmomatic_output_directory}/${Input}${InputFileExtension_rnaSPAdes} ILLUMINACLIP:${Adapters_file}:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25" >> ${Commands_check_log}

         ${trimmomaticPROG} SE -threads ${Num_threads} -trimlog ${Trimmomatic_output_directory}/${Input}.logfile -summary ${Trimmomatic_output_directory}/${Input}_statssummary.txt ${InputPath1}/${Input}.${InputFileExtension_FASTQC_and_Trimmomatic} ${Trimmomatic_output_directory}/${Input}${InputFileExtension_rnaSPAdes} ILLUMINACLIP:${Adapters_file}:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25
             
         mv ${Trimmomatic_output_directory}/${Input}${InputFileExtension_rnaSPAdes} ${Trimmomatic_output_directory}/${ProjectName}/${Input}${InputFileExtension_rnaSPAdes}
         mv ${Trimmomatic_output_directory}/${Input}_statssummary.txt ${Trimmomatic_output_directory}/${ProjectName}/${Input}_statssummary.txt

         rm ${Trimmomatic_output_directory}/${Input}.logfile

         touch ${Trimmomatic_output_directory}/trimmomatic.ok
         echo "Trimmomatic complete for ${Input}" >> ${Pipeline_Logfile}
         ${D_And_T} >> ${Pipeline_Logfile}
      done    

   else
      echo "Script will NOT run Trimmomatic" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi

#----------------------------------------------------------------------#
#                      4. Assembly of contigs                          #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#                      4. Assembly of contigs                          #\n#----------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   # rnaSPAdes
   if [ ${Run_rnaSPAdes} = YES ]
   then 
      echo "Script will run rnaSPAdes" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}

      mkdir ${rnaSPAdes_output_directory}/${ProjectName}

      # To create a string (horizontal list) which has all the sample names with spaces inbetween each name 
      # if you want to do all of this in one line and not have any file remaining afterwards
      InputShortName=($(sed -e "s+${InputFileExtension_rnaSPAdes}++g" <(echo $(find ${InputPath_rnaSPAdes}/${ProjectName} -type f -name "*${InputFileExtension_rnaSPAdes}" -printf "%f "))))
   
      for Input in ${InputShortName[*]} 
      do 
         echo "###Input sample for rnaSPAdes is \"${Input}\"" >> ${Pipeline_Logfile}
         ${D_And_T} >> ${Pipeline_Logfile}
         ##rnaSPAdes
         echo "##rnaSPAdes" >> ${Pipeline_Logfile}
         #rnaSPAdes command
         echo "${SPAdesPROG} --rna -m ${Amount_of_RAM} -t ${Num_threads} -s ${InputPath_rnaSPAdes}/${ProjectName}/${Input}${InputFileExtension_rnaSPAdes} -o ${rnaSPAdes_output_directory}/${ProjectName}/${Input}_rnaSPAdesoutput" >> ${Commands_check_log}

         ${SPAdesPROG} --rna -m ${Amount_of_RAM} -t ${Num_threads} -s ${InputPath_rnaSPAdes}/${ProjectName}/${Input}${InputFileExtension_rnaSPAdes} -o ${rnaSPAdes_output_directory}/${ProjectName}/${Input}_rnaSPAdesoutput

      ${D_And_T} >> ${Pipeline_Logfile}
      done
   else
      echo "Script will NOT run rnaSPAdes" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi

#----------------------------------------------------------------------#
#                  5. Renaming the contig names                        #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#                  5. Renaming the contig names                        #\n#----------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   # Renaming the nodes to include the project name and sample name
   if [ ${Change_contig_names} = YES ]
   then 
      echo "Script will change the assembled contig names using sed" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}

      mkdir ${renamed_assembly_outputs}_relabelled
      mkdir ${renamed_assembly_outputs}_Not_relabelled

      mkdir ${renamed_assembly_outputs}_relabelled/${ProjectName}_samples_${assembler_used}
      labelled=${renamed_assembly_outputs}_relabelled/${ProjectName}_samples_${assembler_used}
         
      mkdir ${renamed_assembly_outputs}_Not_relabelled/${ProjectName}_samples_${assembler_used}
      Not_labelled=${renamed_assembly_outputs}_Not_relabelled/${ProjectName}_samples_${assembler_used}


      # To create a string (horizontal list) which has all the sample names with spaces inbetween each name 
      # if you want to do all of this in one line and not have any file remaining afterwards
      InputShortName=($(sed -e "s+${InputFileExtension_rnaSPAdes}++g" <(echo $(find ${InputPath_rnaSPAdes}/${ProjectName} -type f -name "*${InputFileExtension_rnaSPAdes}" -printf "%f "))))
   
      for Input in ${InputShortName[*]} 
      do 
         echo "###Changing contig names of sample \"${Input}\"" >> ${Pipeline_Logfile}
         ${D_And_T} >> ${Pipeline_Logfile}
         ##sed
         echo "##Changing contig name" >> ${Pipeline_Logfile}
         #Changing name commands
         echo "cp ${rnaSPAdes_output_directory}/${ProjectName}/${Input}_rnaSPAdesoutput/${assembler_fasta_PREFIX}${Diamond_InputFileExtension} ${labelled}/${assembler_fasta_PREFIX}_${Input}${Diamond_InputFileExtension}" >> ${Commands_check_log}

         cp ${rnaSPAdes_output_directory}/${ProjectName}/${Input}_rnaSPAdesoutput/${assembler_fasta_PREFIX}${Diamond_InputFileExtension} ${labelled}/${assembler_fasta_PREFIX}_${Input}${Diamond_InputFileExtension}
         
         sed_renaming_output_type_name=${assembler_fasta_PREFIX}_${Input}${Diamond_InputFileExtension}
         
         echo "sed -ie "s+${contig_name_prefix}+${ProjectName}_${Input}_${contig_name_prefix}+g" ${labelled}/${sed_renaming_output_type_name}" >> ${Commands_check_log}
         sed -ie "s+${contig_name_prefix}+${ProjectName}_${Input}_${contig_name_prefix}+g" ${labelled}/${sed_renaming_output_type_name}
         
         echo "mv ${labelled}/${sed_renaming_output_type_name}e ${Not_labelled}/${sed_renaming_output_type_name}" >> ${Commands_check_log}
         mv ${labelled}/${sed_renaming_output_type_name}e ${Not_labelled}/${sed_renaming_output_type_name}
         
         ${D_And_T} >> ${Pipeline_Logfile}
      done
   else
      echo "Script will NOT change the assembled contig names" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi

#----------------------------------------------------------------------#
#     6. First classification using diamond to obtain viral hits       #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#     6. First classification using diamond to obtain viral hits       #\n#----------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   # Diamond virus screening of assembled contigs
   if [ ${Run_diamond_spades_virus} = YES ]
   then 
      echo "Script will run Diamond BLASTx using virus database" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}
      
      #Making directories
      mkdir ${Diamond_output_directory}/First_classification
      DIAMOND_Assembler_BLAST_Output_Directory=${Diamond_output_directory}/First_classification/${ProjectName}_diamondoutput_${virus_database_shortname}
      mkdir ${DIAMOND_Assembler_BLAST_Output_Directory}
      
      #This variable is called again here incase the renaming contig comand was not run and the variable is not defined beforehand
      labelled=${renamed_assembly_outputs}_relabelled/${ProjectName}_samples_${assembler_used}
      
      # To create a string (horizontal list) which has all the sample names with spaces inbetween each name 
      # if you want to do all of this in one line and not have any file remaining afterwards
      InputShortName=($(sed -e "s+${assembler_fasta_PREFIX}_++g" <(echo $(sed -e "s+${Diamond_InputFileExtension}++g" <(echo $(find ${labelled} -type f -name "*${Diamond_InputFileExtension}" -printf "%f "))))))
   
      for Input in ${InputShortName[*]} 
      do 
         echo "###Input sample for Diamond BLASTx is \"${Input}\"" >> ${Pipeline_Logfile}
         ${D_And_T} >> ${Pipeline_Logfile}
         ##Diamond
         echo "##Diamond" >> ${Pipeline_Logfile}
         #Diamond command
         echo "${diamondProg} blastx --threads ${Num_threads} --range-culling --block-size ${DiamondBlocksize} --max-target-seqs 1 -F 15 --matrix BLOSUM62 -f 100 --salltitles --db ${Virus_DatabaseForDIAMOND} -q ${labelled}/${assembler_fasta_PREFIX}_${Input}${Diamond_InputFileExtension} -o ${DIAMOND_Assembler_BLAST_Output_Directory}/${Input}.daa" >> ${Commands_check_log}

         ${diamondProg} blastx --threads ${Num_threads} --range-culling --block-size ${DiamondBlocksize} --max-target-seqs 1 -F 15 --matrix BLOSUM62 -f 100 --salltitles --db ${Virus_DatabaseForDIAMOND} -q ${labelled}/${assembler_fasta_PREFIX}_${Input}${Diamond_InputFileExtension} -o ${DIAMOND_Assembler_BLAST_Output_Directory}/${Input}.daa
         
         ${D_And_T} >> ${Pipeline_Logfile}
                  
         echo "Importing DIAMOND output to readable format" >> ${Commands_check_log}
         echo "${diamondProg} view --threads ${Num_threads} --daa ${DIAMOND_Assembler_BLAST_Output_Directory}/${Input}.daa --out ${DIAMOND_Assembler_BLAST_Output_Directory}/${Input}_DIAMONDOutput.txt" >> ${Commands_check_log}
        
         ${diamondProg} view --threads ${Num_threads} --daa ${DIAMOND_Assembler_BLAST_Output_Directory}/${Input}.daa --out ${DIAMOND_Assembler_BLAST_Output_Directory}/${Input}_DIAMONDOutput.txt
         
         ${D_And_T} >> ${Pipeline_Logfile}
      done
   else
      echo "Script will NOT run Diamond" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi

#----------------------------------------------------------------------#
#            7. Extraction of contigs with viral hits                  #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#            7. Extraction of contigs with viral hits                  #\n#----------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   # Extracting contigs
   if [ ${Extract_assembled_contigs_with_viral_hits} = YES ]
   then 
      echo "Script will run awk, cut and seqkit to extract contigs from assemblies with viral hits" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}
      
      #This variable is called again here incase diamond was not run and the variable is not defined beforehand
      DIAMOND_Assembler_BLAST_Output_Directory=${Diamond_output_directory}/First_classification/${ProjectName}_diamondoutput_${virus_database_shortname}
      
      #This variable is called again here incase the renaming contig comand was not run and the variable is not defined beforehand
      labelled=${renamed_assembly_outputs}_relabelled/${ProjectName}_samples_${assembler_used}
      
      # Making directory to store lists of sequences to extract
      SeqID_for_extractions_directory=${Extracted_output_contigs}/Sequence_lists_for_extractions
      mkdir ${SeqID_for_extractions_directory}
      mkdir ${SeqID_for_extractions_directory}/first_diamond_extraction
      mkdir ${SeqID_for_extractions_directory}/first_diamond_extraction/${ProjectName}

      First_extractions=${SeqID_for_extractions_directory}/first_diamond_extraction/${ProjectName}
      
      # To create a string (horizontal list) which has all the sample names with spaces inbetween each name 
      # if you want to do all of this in one line and not have any file remaining afterwards
      InputShortName=($(sed -e "s+${assembler_fasta_PREFIX}_++g" <(echo $(sed -e "s+${Diamond_InputFileExtension}++g" <(echo $(find ${labelled} -type f -name "*${Diamond_InputFileExtension}" -printf "%f "))))))
   
      for Input in ${InputShortName[*]} 
      do 
         echo "###Input sample for viral hit extraction is \"${Input}\"" >> ${Pipeline_Logfile}
         ${D_And_T} >> ${Pipeline_Logfile}
         ##awk and cut extraction
         echo "## awk and cut extraction" >> ${Pipeline_Logfile}
         # awk and cut extraction command
         echo "awk 'NR==FNR {accession_file[$1]++; next} $2 in accession_file {print $1}' ${Taxa_accessions_list} ${DIAMOND_Assembler_BLAST_Output_Directory}/${Input}_DIAMONDOutput.txt > ${First_extractions}/${Input}_first_diamond_extraction_contig_names.txt" >> ${Commands_check_log} 

         awk 'NR==FNR {accession_file[$1]++; next} $2 in accession_file {print $1}' ${Taxa_accessions_list} ${DIAMOND_Assembler_BLAST_Output_Directory}/${Input}_DIAMONDOutput.txt > ${First_extractions}/${Input}_first_diamond_extraction_contig_names.txt

         echo "${seqkitProg} grep -f ${First_extractions}/${Input}_first_diamond_extraction_contig_names.txt ${labelled}/${assembler_fasta_PREFIX}_${Input}${Diamond_InputFileExtension} > ${First_extractions}/${Input}_first_extracted_viral_hit_contigs${Diamond_InputFileExtension}" >> ${Commands_check_log}  
      
        ${seqkitProg} grep -f ${First_extractions}/${Input}_first_diamond_extraction_contig_names.txt ${labelled}/${assembler_fasta_PREFIX}_${Input}${Diamond_InputFileExtension} > ${First_extractions}/${Input}_first_extracted_viral_hit_contigs${Diamond_InputFileExtension}

         ${D_And_T} >> ${Pipeline_Logfile}
      done
   else
      echo "Script will not extract contigs from assemblies with viral hits" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi


#----------------------------------------------------------------------#
#            8. Krona display of contigs with viral hits               #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#            8. Krona display of contigs with viral hits               #\n#----------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   # Krona chart generation for diamond output
   if [ ${Run_Krona_on_assembled_contigs} = YES ]
   then 
      echo "Script will run Krona" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}

      #Make directory for sample Krona chart files 
      mkdir ${Krona_output_directory}/First_classification      
      Krona_Assembler_BLAST_Output_Directory=${Krona_output_directory}/First_classification/${ProjectName}_kronaoutput_${virus_database_shortname}
      mkdir ${Krona_Assembler_BLAST_Output_Directory}

      #This variable is called again here incase the renaming contig comand was not run and the variable is not defined beforehand
      labelled=${renamed_assembly_outputs}_relabelled/${ProjectName}_samples_${assembler_used}

      #This variable is called again here incase diamond was not run and the variable is not defined beforehand
      DIAMOND_Assembler_BLAST_Output_Directory=${Diamond_output_directory}/First_classification/${ProjectName}_diamondoutput_${virus_database_shortname}

      echo "source activate ${Krona_environment}" >> ${Commands_check_log}
      source activate ${Krona_environment}
      echo "environment activated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

      # To create a string (horizontal list) which has all the sample names with spaces inbetween each name 
      # if you want to do all of this in one line and not have any file remaining afterwards
      InputShortName=($(sed -e "s+${assembler_fasta_PREFIX}_++g" <(echo $(sed -e "s+${Diamond_InputFileExtension}++g" <(echo $(find ${labelled} -type f -name "*${Diamond_InputFileExtension}" -printf "%f "))))))
   
      for Input in ${InputShortName[*]} 
      do 
         echo "###Input sample for Krona is \"${Input}\"" >> ${Pipeline_Logfile}
         ${D_And_T} >> ${Pipeline_Logfile}
         ##Krona
         echo "##Krona" >> ${Pipeline_Logfile}
         #Krona commands
        echo "Generating Krona chart" >> ${Commands_check_log}
        echo "${KronaBLASTImportPROG} -f -c ${DIAMOND_Assembler_BLAST_Output_Directory}/${Input}_DIAMONDOutput.txt -o ${Krona_Assembler_BLAST_Output_Directory}/${Input}_Krona.html" >> ${Commands_check_log}
        ${KronaBLASTImportPROG} -f -c ${DIAMOND_Assembler_BLAST_Output_Directory}/${Input}_DIAMONDOutput.txt -o ${Krona_Assembler_BLAST_Output_Directory}/${Input}_Krona.html
        
        ${D_And_T} >> ${Pipeline_Logfile}
      done
      
      echo "conda deactivate" >> ${Commands_check_log}
      conda deactivate
      echo "environment deactivated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   else 
      echo "Krona charts will not be made for Diamond virus outputs" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi 
      
   ## Multicharts per bioproject
        
   if [ ${Run_Krona_project_on_assembled_contigs} = YES ]
   then 
      echo "Script will use Krona to make charts for the project" | tee -a ${Pipeline_Logfile} ${Commands_check_log} 
      
      #Activating Krona environment

      echo "source activate ${Krona_environment}" >> ${Commands_check_log}
      source activate ${Krona_environment}
      echo "environment activated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      
      #This variable is called again here incase diamond was not run and the variable is not defined beforehand
      DIAMOND_Assembler_BLAST_Output_Directory=${Diamond_output_directory}/First_classification/${ProjectName}_diamondoutput_${virus_database_shortname}      

      #Make directory for bioproject combined Krona chart files      
      KronaMulti_Assembler_BLAST_OutputDirectory=${Krona_output_directory}/First_classification/${ProjectName}_kronamultioutput_${virus_database_shortname}      
      mkdir ${KronaMulti_Assembler_BLAST_OutputDirectory}
      
      # To create a string (horizontal list) which has all the DIAMOND output sample names with spaces inbetween each name with their file path/address as a prefix
      echo "find ${DIAMOND_Assembler_BLAST_Output_Directory} -type f -name "*_DIAMONDOutput.txt" -printf "%p " > ${KronaMulti_Assembler_BLAST_OutputDirectory}/DIAMONDoutput_${ProjectName}_project_samples_used.txt" >> ${Commands_check_log}

      find ${DIAMOND_Assembler_BLAST_Output_Directory} -type f -name "*_DIAMONDOutput.txt" -printf "%p " > ${KronaMulti_Assembler_BLAST_OutputDirectory}/DIAMONDoutput_${ProjectName}_project_samples_used.txt

find ${DIAMOND_Assembler_BLAST_Output_Directory} -type f -name "*_DIAMONDOutput.txt" -printf "%f " > ${KronaMulti_Assembler_BLAST_OutputDirectory}/DIAMONDoutput_${ProjectName}_project_samples_used_no_address.txt
      
      echo "##Krona multiple sample input" >> ${Pipeline_Logfile}
      ${D_And_T} >> ${Pipeline_Logfile}
      # To turn the input sample string in to a variable to be used in the script
      KronaInputShortName=`cat ${KronaMulti_Assembler_BLAST_OutputDirectory}/DIAMONDoutput_${ProjectName}_project_samples_used.txt`
      echo "${KronaBLASTImportPROG} -f -c ${KronaInputShortName} -o ${KronaMulti_Assembler_BLAST_OutputDirectory}/${ProjectName}_Krona_Multiple_Samples.html" >> ${Commands_check_log}
      
      ${KronaBLASTImportPROG} -f -c ${KronaInputShortName} -o ${KronaMulti_Assembler_BLAST_OutputDirectory}/${ProjectName}_Krona_Multiple_Samples.html 
      ${D_And_T} >> ${Pipeline_Logfile}
      
      echo "conda deactivate" >> ${Commands_check_log}
      conda deactivate
      echo "environment deactivated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
 
   else 
      echo "##Krona will not generate charts for samples combined for ${ProjectName} project" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi


#-------------------------------------------------------------------------------------#
#     9. Second classification using diamond to detect false positives viral hits     #
#-------------------------------------------------------------------------------------#

echo -e "#-------------------------------------------------------------------------------------#\n#     9. Second classification using diamond to detect false positives viral hits     #\n#-------------------------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   # Diamond virus screening of extracted viral hits
   if [ ${Run_diamond_virus_hits_full_db} = YES ]
   then 
      echo "Script will run Diamond BLASTx using full database" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}

      #Making directories
      mkdir ${Diamond_output_directory}/Second_classification
      DIAMOND_virus_final_BLAST_OutputDirectory=${Diamond_output_directory}/Second_classification/${ProjectName}_diamondoutput_${Full_database_shortname}
      mkdir ${DIAMOND_virus_final_BLAST_OutputDirectory}

      #This variable is called again here incase the renaming contig comand was not run and the variable is not defined beforehand
      labelled=${renamed_assembly_outputs}_relabelled/${ProjectName}_samples_${assembler_used}
      
      #These variables are called again here incase the sequence extraction was not run and the variables were not defined beforehand
      SeqID_for_extractions_directory=${Extracted_output_contigs}/Sequence_lists_for_extractions
      First_extractions=${SeqID_for_extractions_directory}/first_diamond_extraction/${ProjectName}
            
      # To create a string (horizontal list) which has all the sample names with spaces inbetween each name 
      # if you want to do all of this in one line and not have any file remaining afterwards
      InputShortName=($(sed -e "s+${assembler_fasta_PREFIX}_++g" <(echo $(sed -e "s+${Diamond_InputFileExtension}++g" <(echo $(find ${labelled} -type f -name "*${Diamond_InputFileExtension}" -printf "%f "))))))   
      for Input in ${InputShortName[*]} 
      do 
         echo "###Input sample for Diamond BLASTx is \"${Input}\"" >> ${Pipeline_Logfile}
         ${D_And_T} >> ${Pipeline_Logfile}
         ##Diamond
         echo "##Diamond" >> ${Pipeline_Logfile}
         #echo "--very-sensitive included" >> ${Pipeline_Logfile}
         #echo "--very-sensitive included" >> ${Commands_check_log}
         #Diamond command
         echo "${diamondProg} blastx --threads ${Num_threads} --range-culling --block-size ${DiamondBlocksize} --max-target-seqs 1 -F 15 --matrix BLOSUM62 -f 100 --salltitles --db ${Full_DatabaseForDIAMOND} -q ${First_extractions}/${Input}_first_extracted_viral_hit_contigs${Diamond_InputFileExtension} -o ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}.daa" >> ${Commands_check_log}

         ${diamondProg} blastx --threads ${Num_threads} --range-culling --block-size ${DiamondBlocksize} --max-target-seqs 1 -F 15 --matrix BLOSUM62 -f 100 --salltitles --db ${Full_DatabaseForDIAMOND} -q ${First_extractions}/${Input}_first_extracted_viral_hit_contigs${Diamond_InputFileExtension} -o ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}.daa
         
         ${D_And_T} >> ${Pipeline_Logfile}
                  
         echo "Importing DIAMOND output to readable format" >> ${Commands_check_log}
         echo "${diamondProg} view --threads ${Num_threads} --daa ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}.daa --out ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}_DIAMONDOutput.txt" >> ${Commands_check_log}
        
         ${diamondProg} view --threads ${Num_threads} --daa ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}.daa --out ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}_DIAMONDOutput.txt
         
         ${D_And_T} >> ${Pipeline_Logfile}
      done
   else
      echo "Script will NOT run Diamond" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi


#----------------------------------------------------------------------#
#           10. Extraction of viral candidate contigs                  #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#           10. Extraction of viral candidate contigs                  #\n#----------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   # Extracting contigs
   if [ ${Extract_viral_candidate_contigs} = YES ]
   then 
      echo "Script will run awk, cut and seqkit to extract viral candidate contigs" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}
      
      #This variable is called again here incase diamond was not run and the variable is not defined beforehand
      DIAMOND_virus_final_BLAST_OutputDirectory=${Diamond_output_directory}/Second_classification/${ProjectName}_diamondoutput_${Full_database_shortname}
      
      #This variable is called again here incase the renaming contig comand was not run and the variable is not defined beforehand
      labelled=${renamed_assembly_outputs}_relabelled/${ProjectName}_samples_${assembler_used}

      #These variables are called again here incase the first sequence extraction was not run and the variables were not defined beforehand
      First_extractions=${SeqID_for_extractions_directory}/first_diamond_extraction/${ProjectName}
      SeqID_for_extractions_directory=${Extracted_output_contigs}/Sequence_lists_for_extractions
      
      # Making directory to store lists of sequences to extract
      mkdir ${SeqID_for_extractions_directory}/second_diamond_extraction
      mkdir ${SeqID_for_extractions_directory}/second_diamond_extraction/${ProjectName}
      mkdir ${Extracted_output_contigs}/virus_candidate_diamond_table
      mkdir ${Extracted_output_contigs}/virus_candidate_diamond_table/${ProjectName}
      
      Second_extractions=${SeqID_for_extractions_directory}/second_diamond_extraction/${ProjectName}
      Virus_candidate_diamond_table_output_directory=${Extracted_output_contigs}/virus_candidate_diamond_table/${ProjectName}
            
      # To create a string (horizontal list) which has all the sample names with spaces inbetween each name 
      # if you want to do all of this in one line and not have any file remaining afterwards
      InputShortName=($(sed -e "s+${assembler_fasta_PREFIX}_++g" <(echo $(sed -e "s+${Diamond_InputFileExtension}++g" <(echo $(find ${labelled} -type f -name "*${Diamond_InputFileExtension}" -printf "%f "))))))
   
      for Input in ${InputShortName[*]} 
      do 
         echo "###Input sample for viral hit extraction is \"${Input}\"" >> ${Pipeline_Logfile}
         ${D_And_T} >> ${Pipeline_Logfile}
         ##awk virus candidate diamond extraction
         echo "## awk virus candidate diamond table extraction" >> ${Pipeline_Logfile}
         # awk virus candidate diamond extraction command
         echo "awk 'NR==FNR {accession_file[$1]++; next} $2 in accession_file {print $0}' ${Taxa_accessions_list} ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}_DIAMONDOutput.txt > ${Virus_candidate_diamond_table_output_directory}/${Input}_virus_candidate_only_DIAMONDOutput.txt" >> ${Commands_check_log}
         
         awk 'NR==FNR {accession_file[$1]++; next} $2 in accession_file {print $0}' ${Taxa_accessions_list} ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}_DIAMONDOutput.txt > ${Virus_candidate_diamond_table_output_directory}/${Input}_virus_candidate_only_DIAMONDOutput.txt
         
         ##awk and cut extraction
         echo "## awk and cut extraction" >> ${Pipeline_Logfile}
         # awk and cut extraction command         
         echo "awk 'NR==FNR {accession_file[$1]++; next} $2 in accession_file {print $1}' ${Taxa_accessions_list} ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}_DIAMONDOutput.txt > ${Second_extractions}/${Input}_virus_candidate_diamond_extraction_contig_names.txt" >> ${Commands_check_log}  

         awk 'NR==FNR {accession_file[$1]++; next} $2 in accession_file {print $1}' ${Taxa_accessions_list} ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}_DIAMONDOutput.txt > ${Second_extractions}/${Input}_virus_candidate_diamond_extraction_contig_names.txt

         echo "${seqkitProg} grep -f ${Second_extractions}/${Input}_virus_candidate_diamond_extraction_contig_names.txt ${First_extractions}/${Input}_first_extracted_viral_hit_contigs${Diamond_InputFileExtension} > ${Second_extractions}/${Input}_virus_candidate_viral_hit_contigs${Diamond_InputFileExtension}" >> ${Commands_check_log}  
      
        ${seqkitProg} grep -f ${Second_extractions}/${Input}_virus_candidate_diamond_extraction_contig_names.txt ${First_extractions}/${Input}_first_extracted_viral_hit_contigs${Diamond_InputFileExtension} > ${Second_extractions}/${Input}_virus_candidate_viral_hit_contigs${Diamond_InputFileExtension}

         ${D_And_T} >> ${Pipeline_Logfile}
      done
   else
      echo "Script will not extract virus candidate contigs from first virus hits" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi


#----------------------------------------------------------------------#
#           11. Krona display of contigs with viral hits               #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#           11. Krona display of contigs with viral hits               #\n#----------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   # Krona chart generation for second diamond output
   if [ ${Run_Krona_on_viral_candidates} = YES ]
   then 
      echo "Script will run Krona" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}

      #Make directory for sample Krona chart files 
      mkdir ${Krona_output_directory}/Second_classification      
      Krona_virus_final_BLAST_OutputDirectory=${Krona_output_directory}/Second_classification/${ProjectName}_kronaoutput_${Full_database_shortname}
      mkdir ${Krona_virus_final_BLAST_OutputDirectory}

      #This variable is called again here incase diamond was not run and the variable is not defined beforehand
      DIAMOND_virus_final_BLAST_OutputDirectory=${Diamond_output_directory}/Second_classification/${ProjectName}_diamondoutput_${Full_database_shortname}

      #This variable is called again here incase the renaming contig comand was not run and the variable is not defined beforehand
      labelled=${renamed_assembly_outputs}_relabelled/${ProjectName}_samples_${assembler_used}
      
      echo "source activate ${Krona_environment}" >> ${Commands_check_log}
      source activate ${Krona_environment}
      echo "environment activated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

      # To create a string (horizontal list) which has all the sample names with spaces inbetween each name 
      # if you want to do all of this in one line and not have any file remaining afterwards
      InputShortName=($(sed -e "s+${assembler_fasta_PREFIX}_++g" <(echo $(sed -e "s+${Diamond_InputFileExtension}++g" <(echo $(find ${labelled} -type f -name "*${Diamond_InputFileExtension}" -printf "%f "))))))
   
      for Input in ${InputShortName[*]} 
      do 
         echo "###Input sample for Krona is \"${Input}\"" >> ${Pipeline_Logfile}
         ${D_And_T} >> ${Pipeline_Logfile}
         ##Krona
         echo "##Krona" >> ${Pipeline_Logfile}
         #Krona commands
        echo "Generating Krona chart" >> ${Commands_check_log}
        echo "${KronaBLASTImportPROG} -f -c ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}_DIAMONDOutput.txt -o ${Krona_virus_final_BLAST_OutputDirectory}/${Input}_Krona.html" >> ${Commands_check_log}
        ${KronaBLASTImportPROG} -f -c ${DIAMOND_virus_final_BLAST_OutputDirectory}/${Input}_DIAMONDOutput.txt -o ${Krona_virus_final_BLAST_OutputDirectory}/${Input}_Krona.html
        
        ${D_And_T} >> ${Pipeline_Logfile}
      done
      
      echo "conda deactivate" >> ${Commands_check_log}
      conda deactivate
      echo "environment deactivated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   else 
      echo "Krona charts will not be made for Diamond virus outputs" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi 
      
   ## Multicharts per bioproject
        
   if [ ${Run_Krona_project_on_viral_candidate_contigs} = YES ]
   then 
      echo "Script will use Krona to make charts for the project" | tee -a ${Pipeline_Logfile} ${Commands_check_log} 
      
      #Activating Krona environment

      echo "source activate ${Krona_environment}" >> ${Commands_check_log}
      source activate ${Krona_environment}
      echo "environment activated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      
      #This variable is called again here incase diamond was not run and the variable is not defined beforehand
      DIAMOND_virus_final_BLAST_OutputDirectory=${Diamond_output_directory}/Second_classification/${ProjectName}_diamondoutput_${Full_database_shortname}
      
      #Make directory for bioproject combined Krona chart files 
      KronaMulti_virus_final_BLAST_OutputDirectory=${Krona_output_directory}/Second_classification/${ProjectName}_kronamultioutput_${Full_database_shortname}
      mkdir ${KronaMulti_virus_final_BLAST_OutputDirectory}
      
      # To create a string (horizontal list) which has all the DIAMOND output sample names with spaces inbetween each name with their file path/address as a prefix
      echo "find ${DIAMOND_virus_final_BLAST_OutputDirectory} -type f -name "*_DIAMONDOutput.txt" -printf "%p " > ${KronaMulti_virus_final_BLAST_OutputDirectory}/DIAMONDoutput_${ProjectName}_project_samples_used.txt" >> ${Commands_check_log}

      find ${DIAMOND_virus_final_BLAST_OutputDirectory} -type f -name "*_DIAMONDOutput.txt" -printf "%p " > ${KronaMulti_virus_final_BLAST_OutputDirectory}/DIAMONDoutput_${ProjectName}_project_samples_used.txt
      
      find ${DIAMOND_virus_final_BLAST_OutputDirectory} -type f -name "*_DIAMONDOutput.txt" -printf "%f " > ${KronaMulti_virus_final_BLAST_OutputDirectory}/DIAMONDoutput_${ProjectName}_project_samples_used_no_address.txt
      
      echo "##Krona multiple sample input" >> ${Pipeline_Logfile}
      ${D_And_T} >> ${Pipeline_Logfile}
      # To turn the input sample string in to a variable to be used in the script
      KronaInputShortName=`cat ${KronaMulti_virus_final_BLAST_OutputDirectory}/DIAMONDoutput_${ProjectName}_project_samples_used.txt`
      
      echo "${KronaBLASTImportPROG} -f -c ${KronaInputShortName} -o ${KronaMulti_virus_final_BLAST_OutputDirectory}/${ProjectName}_Krona_Multiple_Samples.html" >> ${Commands_check_log}
      
      ${KronaBLASTImportPROG} -f -c ${KronaInputShortName} -o ${KronaMulti_virus_final_BLAST_OutputDirectory}/${ProjectName}_Krona_Multiple_Samples.html 

      ${D_And_T} >> ${Pipeline_Logfile}
      
      echo "conda deactivate" >> ${Commands_check_log}
      conda deactivate
      echo "environment deactivated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
 
   else 
      echo "##Krona will not generate charts for samples combined for ${ProjectName} project" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi


#---------------------------------------------------------------------------------#
#           12. Placing all contigs from all samples into one file for CheckV     #
#---------------------------------------------------------------------------------#

echo -e "#---------------------------------------------------------------------------------#\n#           12. Placing all contigs from all samples into one file for CheckV     #\n#---------------------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

   #File concatenation
   if [ ${Do_you_want_to_concatenate_fasta_files_for_CheckV} = YES ]
   then 
      echo "Script will concatenate the viral candidate fasta files for project ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}
         
      #These variables are called again here incase the second sequence extraction was not run and the variables were not defined beforehand
      SeqID_for_extractions_directory=${Extracted_output_contigs}/Sequence_lists_for_extractions
      Second_extractions=${SeqID_for_extractions_directory}/second_diamond_extraction/${ProjectName}
      concatenated_viral_candidates_extractions=${SeqID_for_extractions_directory}/second_diamond_extraction

      
      #Making the directory
      concatenate_OUTPUT=${concatenated_viral_candidates_extractions}/concatenated_viral_can_fasta/${ProjectName}_concatenated_viral_candidates${Diamond_InputFileExtension}
      mkdir ${concatenated_viral_candidates_extractions}/concatenated_viral_can_fasta

      echo "The script will contatenated the input fasta files for easier retrieval of sequences. You can find the output file here: \"${concatenate_OUTPUT}\"" >> ${Pipeline_Logfile}
      
      find ${Second_extractions} -type f -name "*${Diamond_InputFileExtension}" -exec cat {} + > ${concatenate_OUTPUT}
      echo "gzip -c --best ${concatenate_OUTPUT} > ${concatenate_OUTPUT}.gz" >> ${Commands_check_log}
      gzip -c --best ${concatenate_OUTPUT} > ${concatenate_OUTPUT}.gz
      ${D_And_T} >> ${Pipeline_Logfile}
   else
      echo "The script will NOT concatenate the viral candidate fasta files for CheckV" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi

#----------------------------------------------------------------------#
#           13. CheckV assessment of viral contigs candidate           #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#           13. CheckV assessment of viral contigs candidate           #\n#----------------------------------------------------------------------#\n Project is ${ProjectName}" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

  # CheckV
   if [ ${Run_CheckV} = YES ]
   then 
      echo "Script will run CheckV" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
      ${D_And_T} >> ${Pipeline_Logfile}

      #These variables are called again here incase the viral candidate fasta concatenation was not run and the variables were not defined beforehand
      SeqID_for_extractions_directory=${Extracted_output_contigs}/Sequence_lists_for_extractions
      concatenated_viral_candidates_extractions=${SeqID_for_extractions_directory}/second_diamond_extraction
      concatenate_OUTPUT=${concatenated_viral_candidates_extractions}/concatenated_viral_can_fasta/${ProjectName}_concatenated_viral_candidates${Diamond_InputFileExtension}
      

      #Activating environment      
      echo "source activate ${CheckV_environment}" >> ${Commands_check_log}
      source activate ${CheckV_environment}
      echo "environment activated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
          
      echo "###Input sample for CheckV is \"${ProjectName}_concatenated_viral_candidates${Diamond_InputFileExtension}\"" >> ${Pipeline_Logfile}
      ##CheckV
      echo "##CheckV" >> ${Pipeline_Logfile}
      #CheckV command
      echo "${CheckV_PROG} end_to_end ${concatenate_OUTPUT} ${CheckV_output_directory}/${ProjectName}_CheckVout -t ${Num_threads}" >> ${Commands_check_log}
      

      #obtaining number of threads for CheckV to use to avoid over parallelization of the number of threads selected per number of sequences
      
      #num_seq=$(grep ">" ${concatenate_OUTPUT} | wc -l)
      #if [ ${num_seq} -le ${Num_threads} ]
      #then
      #   Num_threads_requested=${num_seq}
      #else
      #   Num_threads_requested=${Num_threads}
      #fi

      ${CheckV_PROG} end_to_end ${concatenate_OUTPUT} ${CheckV_output_directory}/${ProjectName}_CheckVout -d ${CheckV_DB} -t 1
      ${D_And_T} >> ${Pipeline_Logfile}
      
      echo "conda deactivate" >> ${Commands_check_log}
      conda deactivate
      echo "environment deactivated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   else
      echo "Script will NOT run CheckV" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   fi
done


#----------------------------------------------------------------------#
#                14. Krona display of all samples                      #
#----------------------------------------------------------------------#

echo -e "#----------------------------------------------------------------------#\n#                14. Krona display of all samples                      #\n#----------------------------------------------------------------------#" | tee -a ${Pipeline_Logfile} ${Commands_check_log}

# Krona chart for all projects
        
if [ ${Run_Krona_project_on_viral_candidate_contigs_of_all_projects} = YES ]
then 
   echo "Script will use Krona to make charts for all of the projects combined" | tee -a ${Pipeline_Logfile} ${Commands_check_log}  
   ${D_And_T} >> ${Pipeline_Logfile}
   
   #Activating Krona environment

   echo "source activate ${Krona_environment}" >> ${Commands_check_log} 
   source activate ${Krona_environment}
   echo "environment activated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
   
   #This variable is called again here incase the viral candidate contig extraction was not run and the variable is not defined beforehand
   Virus_candidate_diamond_table_output_directory=${Extracted_output_contigs}/virus_candidate_diamond_table

   #Make directory for bioproject combined Krona chart files 
   Krona_all_projects_virus_final_BLAST_Output_Directory=${Krona_output_directory}/Second_classification/${Project_collection}_kronamulti_output_virome
   mkdir ${Krona_all_projects_virus_final_BLAST_Output_Directory}

   echo "##List of samples used printed to ${Krona_all_projects_virus_final_BLAST_Output_Directory}/All_samples_used_for_krona_chart.txt" >> ${Pipeline_Logfile}
   echo "find ${Virus_candidate_diamond_table_output_directory} -type f -name "*_virus_candidate_only_DIAMONDOutput.txt" -printf "%f\n" > ${Krona_all_projects_virus_final_BLAST_Output_Directory}/All_samples_used_for_krona_chart.txt" >> ${Commands_check_log}  
      
   
   find ${Virus_candidate_diamond_table_output_directory} -type f -name "*_virus_candidate_only_DIAMONDOutput.txt" -printf "%f\n" > ${Krona_all_projects_virus_final_BLAST_Output_Directory}/All_samples_used_for_krona_chart.txt
   sed -ie "s+_virus_candidate_only_DIAMONDOutput.txt}++g" ${Krona_all_projects_virus_final_BLAST_Output_Directory}/All_samples_used_for_krona_chart.txt
   rm ${Krona_all_projects_virus_final_BLAST_Output_Directory}/All_samples_used_for_krona_chart.txte

   # To turn the input sample string in to a variable to be used in the script
   ## tell find to print full file path and place this string in the variable, this will get all the sample diamond output file names to be used for krona and 
   ### include the path to these files, which will let the import command find them all
   KronaAllInputShortNames=$(echo $(find ${Virus_candidate_diamond_table_output_directory} -type f -name "*_virus_candidate_only_DIAMONDOutput.txt" -printf "%p "))
   
   echo "##Krona" >> ${Pipeline_Logfile}
   #Krona command 

   ${KronaBLASTImportPROG} -f -c ${KronaAllInputShortNames} -o ${Krona_all_projects_virus_final_BLAST_Output_Directory}/${Project_collection}_Krona_Multiple_Samples.html

   ${D_And_T} >> ${Pipeline_Logfile}
      
   echo "conda deactivate" >> ${Commands_check_log}
   conda deactivate
   echo "environment deactivated" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
else 
   echo "##Krona will not generate charts for all projects combined" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
fi

##################################################################################################################

${D_And_T} | tee -a ${Pipeline_Logfile} ${Commands_check_log}
echo "Pipeline is finished, please check the ${Pipeline_Logfile} and ${Commands_check_log} if you suspect any issues" | tee -a ${Pipeline_Logfile} ${Commands_check_log}
