#!/bin/Rscript
###################################################################################################
# syntax: Rscript techVal_heatmap.R <data.table>
# eg Rscript techVal_heatmap.R merge_calls/delly.tumorSpecific.highConf.cutoff.binaryTree.txt
# Joachim Weischenfeldt, 2013
#
# dataTable format
#   clique DKFZ CNAG OICR SANGER RIKEN
# 1      1    1    1    1      1     1
# 2      2    0    1    0      1     1
# 3      3    0    1    0      1     1
# 4      4    1    1    1      1     1
# 5      5    1    1    1      1     1
# 6      6    1    0    0      1     0
###################################################################################################

require(ggplot2)
require(reshape2)
require(gridExtra)
require(ggdendro)

args <- commandArgs(TRUE)
cat(args)
dataTable <-read.table(args[1], header=TRUE)

head(dataTable)
count = NA
nrow(dataTable)
if ("count" %in% dataTable$clique) {
  count = dataTable[dataTable$clique=="count",]
  dataTable = dataTable[1:(nrow(dataTable)-1),]
}

embl_broad_c = nrow(dataTable[dataTable$embl==1 & dataTable$broad==1 ,])
embl_sanger_c = nrow(dataTable[dataTable$embl==1 & dataTable$sanger==1 ,])
broad_sanger_c = nrow(dataTable[dataTable$broad==1 & dataTable$sanger==1 ,])
embl_broad_sanger_c = nrow(dataTable[dataTable$embl==1 & dataTable$broad==1 & dataTable$sanger==1 ,])

df.m <- melt(dataTable, id.vars=names(dataTable)[1])
summary(df.m)

# plotting

dirName <- sub("(.*).txt","\\1",args, perl=T)
dirName
fileName <- sub(".*/(.*)binaryTree.txt","\\1",args, perl=T)

##
datbin = dataTable[,c(2:ncol(dataTable))]
d = dist(t(datbin), method = "binary")
hc = hclust(d, method="ward")
p1 = ggdendrogram(hc, rotate = FALSE) + labs(title=paste("PCAWG SV caller clustering\n",fileName), size=9) + theme(title = element_text( size=9))

cat (paste( "########### Generating Heatmap for ", fileName," ###########", sep=""))

p2 <- (ggplot(df.m, aes(variable, clique)) + geom_tile(aes(fill = value), colour = "white") + 
        scale_fill_gradient(low = "white", high = "darkseagreen4")) +
 ylim(rev(levels(df.m$clique)))

p2 = p2  +  theme(title = element_text( size=9),
      axis.ticks = element_blank(), 
      axis.text.y =  element_blank(),
      axis.title.x = element_text(colour="#990000", size=14),
      axis.title.y = element_text(colour="#990000", size=14),
      axis.text.x  = element_text(angle=45, vjust=0.5, size=12)) + 
  theme(legend.position="none") +
  labs(title=paste("PCAWG SV concordance\n",
                   "\n\nSVs:, \nBRASS = ", count$brass,  "\nDELLY = ", count$delly, "\n dRANGER = ", count$dranger,
                   "\nSNOWMAN = ", count$snowman,",\nTotal shared SVs = ", nrow(dataTable), sep=""), 
  x="Centers", y="SVs (clusters)" )  

plotFile = paste(dirName, ".heatmap.pdf", sep="")
pdf(file=plotFile, paper='a4r')
grid.arrange(p1, p2, ncol=2)
dev.off()

cat (paste("Output: ", plotFile, "\n\n", sep=""))

