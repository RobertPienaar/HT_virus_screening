#### clear env - load wd ####

rm(list=ls())

setwd("/path/to/directory/housing/gff_files")

{

list_gff <- list.files(pattern="gff$", full.names = TRUE, recursive = TRUE)
#setwd("C:/Users/Utilisateur/Desktop/Insect_doctors_course2/SIP_2022_conference/results/July23_25th2022/Annotated_viruses") # set you won path


annot_name = paste("annot_",gsub("(./)(.*)(.gff)","\\2",list_gff),sep="")

list_annot = list()
for(files in 1:length(list_gff)){
 list_annot[[files]] <- read.delim2(list_gff[files], header =F , sep="\t")
}

names(list_annot) = annot_name

for(df in 1:length(annot_name)){
  list_annot[[df]] <- list_annot[[df]][-c(1:3),]
  colnames(list_annot[[df]]) <- c("contig_name","source","type","start","stop","score","direction",
                            "strand")
  list_annot[[df]] = na.omit(list_annot[[df]])
  
  for(contigs in 1:length(unique(list_annot[[df]][["contig_name"]]))){
    list_annot[[df]][["contig_name"]][list_annot[[df]][["contig_name"]]==unique(list_annot[[df]][["contig_name"]])[contigs]] = paste("contig",contigs,sep=" ")
  }
  
}

#to add colours

T_all_list <- lapply(1:length(list_annot), 
             function(p) {
               T_all <- list_annot[[p]][["type"]]
             }
)

T_uniq <- unique(unlist(T_all_list))


#T_uniq <- unique(types_f(list_annot))
colors <- hcl.colors(length(T_uniq), palette="Zissou 1", alpha=NULL)
col_type <- as.data.frame(cbind(T_uniq,colors))


for(df in 1:length(list_annot)){
list_annot[[df]][["colours"]] <- paste(col_type$colors)[match(list_annot[[df]][["type"]], col_type$T_uniq)]
}

## load your file containing your gene information

#annot <- read.delim2("HiCV1.gff", header =F , sep="\t")
#annot <- annot[-c(1:3),]
#colnames(annot) <- c("contig_name","source","type","start","stop","score","direction",
 #                    "strand","info")
#annot = na.omit(annot)
#for(contigs in 1:length(unique(annot$contig_name))){
 # annot$contig_name[annot$contig_name==unique(annot$contig_name)[contigs]] = paste("contig",contigs,sep=" ")
#}

#annot$contig_name = "HiCV"  # Only use this if there is just one contig
#View(annot)

# if you don't have the length, you can calculate it by max - min + 1


#### plot virus  ####

## first you have to calculate the size of the longest sequence that you want to plot

for(l in 1:length(list_annot)){
  
annot = list_annot[[l]]

max_length = max(annot$stop) - min(annot$start) + 1

# if you have several contigs to plot :
length_contig = c()
for(i in 1:length(unique(annot$contig_name))){ 
  length_contig[i] = max(annot$stop[annot$contig_name==unique(annot$contig_name)[i]]) -
    min(annot$start[annot$contig_name==unique(annot$contig_name)[i]]) + 1
}
max_length_contig = max(length_contig)

annot = annot[which(!annot$type==c("extracted region","Region")),]

## new plot
#tiff(paste(annot_name[l],".tif",sep=""),  units="in", width=20, height=13, res=400, pointsize = 15)
svg(paste(annot_name[l],".svg",sep=""), width=20, height=13, pointsize = 15)


par(mar=c(0,4,1,0.5))
plot.new()

# here you can generates the heights of your sequences for the plot (where it will be plotted)
# I have several sequences/ contigs to plot, so I'm creating a string of heigth values

heights = rev(seq(0.5*(1/(log2(length(unique(annot$contig_name)))+1)), by = 0.3*(3/length(unique(annot$contig_name))),
              length = length(length_contig))) # to change = contig height
# *(1/length(unique(annot$contig_name)))
# here is to build your contig segments on the plot

for(j in 1:length(length_contig)){
  segments(0, heights[j], length_contig[j]/max_length_contig, heights[j], lwd=3)
}

# I'm building a ladder to indicate the sizes 
if(max_length_contig > 10000){b=5000}else if(max_length_contig <= 10000){b=500}

segments(0,1, max_length_contig/max_length_contig,1, col="black", lwd=2)
ladder_range=seq(0, by = b, max_length_contig)
text(x=ladder_range/max_length_contig, y=1.03, labels=format(ladder_range, big.mark = ","), col="black", cex=1, font=2) # to change = ladder numbers
text(x=ladder_range/max_length_contig, y=1.01, labels="-", col="black", cex=1, srt=90)


# now, I want to name my sequences

cluster_names = unique(annot$contig_name) 
for(k in 1:length(cluster_names)){
  mtext(text=strsplit(sub("annot_",replacement="",annot_name), "_")[[l]][1], side=2, at=heights[k]+0.03, cex=1.5, line=-1.5, outer=F,las=2)
  mtext(text=cluster_names[k], side=2, at=heights[k], cex=1.5, line=-1.5, outer=F,las=2) # to change = contig name
  mtext(text=paste("(",format(length_contig[k], big.mark = ",")," bp)", sep=""), 
        side=2, at=heights[k]-0.03, cex=1.1, line=-1.5, outer=F,las=2) # to change = contig size
}


# We need to define the length of the triangle portion of each arrow and the height of our gene shapes:
arrowLen = 0.01;
boxHeight = 0.03
#boxHeight = 0.04*(2/length(unique(annot$contig_name)));

# here the minimal gene positions found for each clusters 

min_contig = c()
for(i in 1:length(unique(annot$contig_name))){ 
  min_contig[i] = min(annot$start[annot$contig_name==unique(annot$contig_name)[i]])
}

# define nested sequences
annot$size <- annot$stop - annot$start
annot$nest <-c(rep(c("yes"), length(annot$size)))

for(a in 1:length(unique(annot$contig_name))){
  for(c in 1:length(annot$size[annot$contig_name==unique(annot$contig_name)[a]])){
  if((annot$size[annot$contig_name==unique(annot$contig_name)[a]][c] >= max(annot$size[annot$contig_name==unique(annot$contig_name)[a]]))=="TRUE"){
    annot$nest[annot$contig_name==unique(annot$contig_name)[a]][c] <- c("no")
  }
}
}
  
## here is the big loop to create the shapes
# for each contig
for(j in 1:length(unique(annot$contig_name))){
  # we calculate the coordinates of the shapes for each genes  
  for (i in 1:length(annot$start[annot$contig_name == unique(annot$contig_name)[j]])){
    left <- (annot$start[annot$contig_name == unique(annot$contig_name)[j]][i]-min_contig[j])/max_length_contig
    right <- (annot$stop[annot$contig_name == unique(annot$contig_name)[j]][i]-min_contig[j])/max_length_contig
    if (annot$nest[annot$contig_name == unique(annot$contig_name)[j]][i]=="no"){
    height <- heights[j]
    print(heights[j])
    }else if(annot$nest[annot$contig_name == unique(annot$contig_name)[j]][i]=="yes"){
      height <- heights[j]-0.045
      print( heights[j]-0.045)
    }
# generating a color string to fill the shapes     
        geneFillColor=as.character(annot$colours[annot$contig_name == unique(annot$contig_name)[j]][i]);
      
# here we draw the shapes according to their position and the direction of the genes
# strand 1 (forward)    
    if (annot$direction[annot$contig_name == unique(annot$contig_name)[j]][i] == "+"){
      arrowStart = max (left, right - arrowLen)
      polygon (c(left, arrowStart, right, arrowStart, left), c(height-boxHeight/2, height-boxHeight/2,
                                                               height, height+boxHeight/2, height+boxHeight/2), col = geneFillColor, lwd = 2)
# strand -1 (reverse)  
    }else if (annot$direction[annot$contig_name == unique(annot$contig_name)[j]][i] == "-"){
      arrowStart = min (right, left + arrowLen)
      polygon (c(left, arrowStart, right, right, arrowStart), c(height, height-boxHeight/2,
                                                                height-boxHeight/2, height+boxHeight/2, height+boxHeight/2), col = geneFillColor, lwd = 2)
    }

  }
} 


dev.off()


system2("open", args = paste(annot_name[l],".svg",sep=""))


}

}

