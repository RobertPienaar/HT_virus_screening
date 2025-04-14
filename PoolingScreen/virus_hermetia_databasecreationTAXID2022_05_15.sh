#!/bin/bash

eval "$(conda shell.bash hook)"

date
## Script for database creation

#Based on your preferred TaxaIDs this script will create a mini BLAST database by extracting from either the NT or NR non-preformatted fasta databases.
# Tip: remember whether you upload the script or use nano or vim to create the file and paste the text in. Always use the command "chmod a+x filename.sh" to make the script executable.
# Tip: This script was written for bash (#!/bin/bash) and not for shell (#!/bin/sh), so be mindful of the syntax if you want to change it.

############## Variables to change ##############

# To call the correct version of taxonkit and seqkit that you have installed, place the file path or your command used to run the software from anywhere on your system.
taxonkitPROG=taxonkit
seqkitPROG=seqkit

## File paths
# Two were place incase the user had limited space and wanted to store the temporary database files in a separate directory. If it doesn't matter, then keep both directories the same for easy use of the script.

taxonPATH=/home/Dataanalysis/databases/taxonkit/downloaded/generated_dbs
NCBIdatabasePATH=/home/Dataanalysis/databases/taxonkit/downloaded

Do_you_want_to_move_output_files_to_taxonPATH=YES

## Tax ID number
# Tip: If you want to use multiple TaxIDs, then you can separate them with a comma in the variable below (e.g. taxaIDnumber=2,4751,10239,50557,2157)

taxaIDnumber=10239,343581

Generate_taxID_file_using_taxonkit=YES


## File names
# Output file from Taxonkit extraction

taxIDfile=virus_hermetia_taxidlist_2022_05_23-taxid.txt

# Accession database used  

accessiondb=prot.accession2taxid.FULL_2022_05_15.gz

# NCBI protein or nucleotide fasta file (nr.gz or nt.gz)

UnfilteredDBgz=original/nr_2022_05_15.gz
#UnfilteredDBextracted=nr_2022_05_15.fa

# Unprocessed extracted accession output file

accession2taxidfile=nr_virus_hermetia_id_acc2taxid.txt

# Processed extracted acession output file

accessionFILE=nr_virus_hermetia_id.acc.txt

# Final database to use for analyses

FilteredDBfasta=nr_virus_hermetia_db_2022_05_15.fa.gz


## Nucleotide or Protein database
# For "To_download_database", indicate either "YES" or "NO" depending on in you want the script to actually download to predetermined databases.
# Tip: You can choose here if you want to download the databases for protein (nr) or for the nucleotide (nt) sequences. Simply type "nt" or "nr" into the variable like so "Downloading_nt_or_nr_database=nr"

To_download_database=NO
Downloading_nt_or_nr_database=nr

############## End of variables to change ##############

#cd ${NCBIdatabasePATH}/

###### fetching TaxIDs 

if [ ${Generate_taxID_file_using_taxonkit} = YES ]
then 
   echo "Creating taxID file ${taxIDfile} using taxID/s ${taxaIDnumber}"
   date 
   ${taxonkitPROG} list --ids ${taxaIDnumber} --indent "" > ${NCBIdatabasePATH}/${taxIDfile}
else 
   echo "Using custom taxID file"
fi 

###### check taxID file##
wc -l ${NCBIdatabasePATH}/${taxIDfile}

###### Downloading nt or nr accession2taxid file and extraction accession numbers

if [ ${To_download_database} = YES ]
then 
   echo "Script will download the nr or nt database"
   date 
   if [ ! ${Downloading_nt_or_nr_database} = nr ]
   then 
      echo "Downloading nucl_gb.accession2taxid.gz"
      wget --output-document ${NCBIdatabasePATH}/${accessiondb} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
   else 
      echo "Downloading prot.accession2taxid.FULL.gz"
      wget --output-document ${NCBIdatabasePATH}/${accessiondb} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
   fi 

  ###### Downloading nt or nr database fasta file
   date
   if [ ! ${Downloading_nt_or_nr_database} = nr ]
   then 
      echo "Downloading nt.gz"
      wget --output-document ${NCBIdatabasePATH}/${UnfilteredDBgz} https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
   else 
      echo "Downloading nr.gz"
      wget --output-document ${NCBIdatabasePATH}/${UnfilteredDBgz} https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
   fi 
else
   echo "Script will NOT download the nr or nt database, please ensure that you have indicated where to find your predownloaded database in the NCBIdatabasePATH variable and put the name with the file extension in the UnfilteredDBgz variable and make sure that it is in a .gz archive"
fi

###### obtaining accessions using taxIDs 

gunzip -c ${NCBIdatabasePATH}/${accessiondb} | csvtk grep -t -f taxid -P ${NCBIdatabasePATH}/${taxIDfile} | csvtk -t cut -f accession.version,taxid | sed 1d > ${NCBIdatabasePATH}/${accession2taxidfile}
cut -f 1 ${NCBIdatabasePATH}/${accession2taxidfile} > ${NCBIdatabasePATH}/${accessionFILE}

###### Check accession file
wc -l ${NCBIdatabasePATH}/${accessionFILE}

###### Filter reads from nt according to ur text file created earlier##
gzip -c -d ${NCBIdatabasePATH}/${UnfilteredDBgz} | ${seqkitPROG} grep -f ${NCBIdatabasePATH}/${accessionFILE} - > ${NCBIdatabasePATH}/${FilteredDBfasta}

echo "Database generated for taxaID/s ${taxaIDnumber}"
date

if [ ${Do_you_want_to_move_output_files_to_taxonPATH} = YES ]
then 
   echo "Moving files to preferred directory ${taxonPATH}"
   mv ${NCBIdatabasePATH}/${accessionFILE} ${taxonPATH}/${accessionFILE}.txt
   mv ${NCBIdatabasePATH}/${taxIDfile} ${taxonPATH}/${taxIDfile}.txt
   mv ${NCBIdatabasePATH}/${accession2taxidfile} ${taxonPATH}/${accession2taxidfile}.txt
   mv ${NCBIdatabasePATH}/${FilteredDBfasta} ${taxonPATH}/${FilteredDBfasta}
else 
   echo "output files will remain in ${NCBIdatabasePATH}"
fi
date
echo "Script complete"
