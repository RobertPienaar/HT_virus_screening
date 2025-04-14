#!/bin/bash

eval "$(conda shell.bash hook)"

# Trimming adapters and performing a QC on illumina datasets
# Written by Robert D. Pienaar
# Version: 3
# Last updated: 2024-12-06 
# Made usable with just one conda environment


###### __________ Importing config file and testing it before proceeding _______ ######

source ./virus_mapping_screen_config_file.yaml

Config_check=$(< <(grep -v "#" ./virus_mapping_screen_config_check.txt | cut -f1))

for VAR_NAME in ${Config_check[*]}
do 
	VAR_NAME_print=${VAR_NAME}
	VAR_NAME=${!VAR_NAME}
	if [[ ! -z ${VAR_NAME} ]]
	then
		echo "Variable: ${VAR_NAME_print}"
		echo "What variable is printing: ${VAR_NAME}"
		echo "-----------------"
		echo "Status: Object inside"
		echo "_________________"
	else
		echo "Variable: ${VAR_NAME_print}"
		echo "What variable is printing: ${VAR_NAME}"
		echo "-----------------"
		echo "Status: EMPTY"
		echo "_________________"
		exit
	fi
done

# Formatting date and time printing in script and file names
D_And_T="date --universal +%A__%d-%B-%Y__%H:%M:%S"
Date_of_pipeline_run="date +%d_%B_%Y"

Pipeline_Logfile=${Pipeline_Logfile_prefix}_`${Date_of_pipeline_run}`.log

echo "${Pipeline_Logfile}"

####################################################################

# Activating conda environment if needed
if [ ${Use_CONDA} = YES ]
then
#Activating mapper environment
source activate ${mapping_env_dir}
else
   echo "Conda environments not activated"
fi

${D_And_T} > ${Pipeline_Logfile}

# Compressing fastq samples
for ProjectName in ${BioprojectName[*]} 
do 
   InputPath_PROJ=${InputPath}/${ProjectName}
   InputShortName=($(sed -e "s+${Input_File_Extension}++g" <(echo $(find ${InputPath_PROJ} -type f -name "*${Input_File_Extension}" -printf "%f "))))
   
   for Input in ${InputShortName[*]}
   do
      ##             Archiving .gz
      echo "##               Archiving using pigz" >> ${Pipeline_Logfile}
      echo "###      Input samples for archiving are \"${Input}\"" >> ${Pipeline_Logfile}
      ${D_And_T} >> ${Pipeline_Logfile}
	  
      if [ ${REACTION} = PAIRED ]
      then 
         ${compression_PROG} -9 -k -p${Num_threads} ${InputPath_PROJ}/${Input}_1.fastq
         ${compression_PROG} -9 -k -p${Num_threads} ${InputPath_PROJ}/${Input}_2.fastq
      else
         ${compression_PROG} -9 -k -p${Num_threads} ${InputPath_PROJ}/${Input}.fastq
      fi
	  
      ${D_And_T} >> ${Pipeline_Logfile}
      echo "# \"${Input}\" archived with pigz" >> ${Pipeline_Logfile}
   done
   echo "##          File tree for \"${InputPath}\"" >> ${Pipeline_Logfile}
   echo "`tree -sh`" >> ${Pipeline_Logfile}
   echo "#_______________________________________________________________________________#" >> ${Pipeline_Logfile}
done

# fastp
mkdir ${trim_OUTPATH}

for ProjectName in ${BioprojectName[*]} 
do 
      mkdir ${trim_OUTPATH}/${ProjectName}
      
	  archive_ext=${Input_File_Extension}.gz
      trimmed_out=${trim_OUTPATH}/${ProjectName}
      InputPath_PROJ=${InputPath}/${ProjectName}
      InputShortName=($(sed -e "s+${archive_ext}++g" <(echo $(find ${InputPath_PROJ} -type f -name "*${archive_ext}" -printf "%f "))))

      for Input in ${InputShortName[*]}
      do
         echo "###      Input samples for fastp are \"${Input}\"" >> ${Pipeline_Logfile}
         ${D_And_T} >> ${Pipeline_Logfile}

         if [ ${REACTION} = PAIRED ]
         then 
            # Paired-end
            ${Trimming_PROG} --thread ${Num_threads} -i ${InputPath_PROJ}/${Input}_1.fastq.gz -I ${InputPath_PROJ}/${Input}_2.fastq.gz -o ${trimmed_out}/${Input}_trimmed_paired1.fastq.gz -O ${trimmed_out}/${Input}_trimmed_paired2.fastq.gz --unpaired1 ${trimmed_out}/${Input}_trimmed_unpaired1.fastq.gz --unpaired2 ${trimmed_out}/${Input}_trimmed_unpaired2.fastq.gz ${fastp_param_PE} --html ${trimmed_out}/${Input}_report.html
         else
            # Single-end
            ${Trimming_PROG} --thread ${Num_threads} -i ${InputPath_PROJ}/${Input}.fastq.gz -o ${trimmed_out}/${Input}_trimmed.fastq.gz ${fastp_param_SE} --html ${trimmed_out}/${Input}_report.html
         fi
         ${D_And_T} >> ${Pipeline_Logfile}
      done
done

#----------------------------------------------------------------------------------------------------------------------#
# Mapping script

echo "This pipeline mapped samples to viruses" > ${Pipeline_Logfile}
${D_And_T} >> ${Pipeline_Logfile}

###################################
#                                 #
# 1. Building indexes for mapping #
#                                 #
###################################

if [ ${Build_indexes_for_mapping} = YES ]
then 
   for vir_ref in ${virus_reference_sequence[*]}
   do 
     # Generating mapper index
     ${mapping_PROG}-build ${input_reference_dir}/${vir_ref}.fasta ${mapper_index_dir}/${vir_ref}

     # Generating fai index 
     ${samtools_PROG} faidx ${input_reference_dir}/${vir_ref}.fasta
   done
else
   echo "No new indexes generated for mapping" >> ${Pipeline_Logfile}
fi

##############################################################################################################
#                                                                                                            #
# 2. creating project directories for mapped outputs for each virus reference sequence according to projects #
#                                                                                                            #
##############################################################################################################

mkdir ${mapping_outputs_dir}

for vir_ref in ${virus_reference_sequence[*]}
do
  mkdir ${mapping_outputs_dir}/${vir_ref}
  Ref_mapped_out=${mapping_outputs_dir}/${vir_ref}
  for ProjectName in ${BioprojectName[*]}
  do 
    mkdir ${Ref_mapped_out}/${ProjectName}
  done
done

#####################
#                   #
# 3. mapping script #
#                   #
#####################

if [ ${REACTION} = PAIRED ]
then 
   # Paired-end

   # By virus reference sequence
   for vir_ref in ${virus_reference_sequence[*]}
   do 
     # By project
     for ProjectName in ${BioprojectName[*]}
     do 
       # By samples in project
       trimmed_out=${trim_OUTPATH}/${ProjectName}
       InputShortName=($(sed -e "s+_trimmed_paired1.fastq.gz++g" <(echo $(find ${trimmed_out} -type f -name "*_trimmed_paired1.fastq.gz" -printf "%f "))))
       Ref_mapped_out=${mapping_outputs_dir}/${vir_ref}/${ProjectName}
       for Input in ${InputShortName[*]}
       do 
         # Mapping the raw reads 
         ${mapping_PROG} -p ${Map_threads} -q --no-unal --phred33 -x ${mapper_index_dir}/${vir_ref} -1 ${trimmed_out}/${Input}_trimmed_paired1.fastq.gz -2 ${trimmed_out}/${Input}_trimmed_paired2.fastq.gz -S ${Ref_mapped_out}/${Input}_mapped.sam
         # Samtools import
         ${samtools_PROG} import ${Input_reference_dir}/${vir_ref}.fasta.fai ${Ref_mapped_out}/${Input}_mapped.sam ${Ref_mapped_out}/${Input}_mapped.bam
         # Sorting with samtools
         ${samtools_PROG} sort ${Ref_mapped_out}/${Input}_mapped.bam -o ${Ref_mapped_out}/${Input}_mappedsorted.bam
         # Creating index for the bam file using samtools
         ${samtools_PROG} index ${Ref_mapped_out}/${Input}_mappedsorted.bam ${Ref_mapped_out}/${Input}_mappedsorted.bai
         ${samtools_PROG} view -b -F 4 ${Ref_mapped_out}/${Input}_mappedsorted.bam > ${Ref_mapped_out}/${Input}_mappedsortedfiltered.bam
         # Creating index for the bam file using samtools
         ${samtools_PROG} index ${Ref_mapped_out}/${Input}_mappedsortedfiltered.bam ${Ref_mapped_out}/${Input}_mappedsortedfiltered.bai
         # Removing the complete mapped Sam and Bam files to save space
         rm ${Ref_mapped_out}/${Input}_mapped.sam
         rm ${Ref_mapped_out}/${Input}_mapped.bam
       done
     done
   done

else 

   # Single-end

   # By virus reference sequence
   for vir_ref in ${virus_reference_sequence[*]}
   do 
     # By project
     for ProjectName in ${BioprojectName[*]}
     do 
       # By samples in project
       trimmed_out=${trim_OUTPATH}/${ProjectName}
       InputShortName=($(sed -e "s+_trimmed.fastq.gz++g" <(echo $(find ${trimmed_out} -type f -name "*_trimmed.fastq.gz" -printf "%f "))))
       Ref_mapped_out=${mapping_outputs_dir}/${vir_ref}/${ProjectName}
       for Input in ${InputShortName[*]}
       do 
         # Mapping the raw reads 
         ${mapping_PROG} -p ${Map_threads} -q --no-unal --phred33 -x ${mapper_index_dir}/${vir_ref} -U ${trimmed_out}/${Input}_trimmed.fastq.gz -S ${Ref_mapped_out}/${Input}_mapped.sam
         # Samtools import
         ${samtools_PROG} import ${Input_reference_dir}/${vir_ref}.fasta.fai ${Ref_mapped_out}/${Input}_mapped.sam ${Ref_mapped_out}/${Input}_mapped.bam
         # Sorting with samtools
         ${samtools_PROG} sort ${Ref_mapped_out}/${Input}_mapped.bam -o ${Ref_mapped_out}/${Input}_mappedsorted.bam
         # Creating index for the bam file using samtools
         ${samtools_PROG} index ${Ref_mapped_out}/${Input}_mappedsorted.bam ${Ref_mapped_out}/${Input}_mappedsorted.bai
         ${samtools_PROG} view -b -F 4 ${Ref_mapped_out}/${Input}_mappedsorted.bam > ${Ref_mapped_out}/${Input}_mappedsortedfiltered.bam
         # Creating index for the bam file using samtools
         ${samtools_PROG} index ${Ref_mapped_out}/${Input}_mappedsortedfiltered.bam ${Ref_mapped_out}/${Input}_mappedsortedfiltered.bai
         # Removing the complete mapped Sam and Bam files to save space
         rm ${Ref_mapped_out}/${Input}_mapped.sam
         rm ${Ref_mapped_out}/${Input}_mapped.bam
       done
     done
   done
fi

conda deactivate

echo "Script has finished" >> ${Pipeline_Logfile}
${D_And_T} >> ${Pipeline_Logfile}
