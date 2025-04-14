Here are the scripts for the high throughput virus screening by mapping.


The recommended computing resources would be to use a high resource computing cluster. Typlically this pipeline is run using 24 cores and 100 GB of ram. However the size of the dataset and reference sequence will influence this. For instance, the visualisations can easily be run on a laptop with 16GB of ram and 8 cores/threads. The final visualisation of the BAM files in Rstudio on a laptop may take some time (more than 30 minutes) to complete per reference sequence if trying to visualise more than 50 datasets together, depending on the size of the BAM files. 
