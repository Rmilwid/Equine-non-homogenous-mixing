#this script take an xlsx file in which the 1B-1n cells have the headers vn, for vaccination threshold
# and the A1-An cells have the row headers qn. the cells in teh matrix should contain the total prevalence for each model at
#each threshold.
# this code turns that matrix into a heatmap
library(lattice)
library(xlsx)
library(gplots)
library(raster)
library(rasterVis)


fp1="C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/network with degree 4/nVacVSnquar.tiff"
trial=read.xlsx("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/network with degree 4/varying intervention summary matrices.xlsx",sheetIndex = 1)
mat_data<-data.matrix(trial[,2:ncol(trial)])
rnames<-trial[,1]
rnames<-gsub("q","",rnames)
rownames(mat_data)<-rnames
colnames(mat_data)<-gsub("v","",colnames(mat_data))
tiff(fp1,height=180, width=180,units='mm', compression="lzw",res=300)
# 
# heatmap.2(mat_data,
#           notecol="black",      # change font color of cell labels to black
#           density.info="none",  # turns off density plot inside color legend
#           trace="none",         # turns off trace lines inside the heat map
#           col=colorRampPalette(c("white","grey","black"))(100),
#           dendrogram="none",
#           Rowv=F,
#           Colv=F,
#           revC=F,
#           cexCol=1.5,
#           cexRow=1.5,
#           srtCol = 0,
#           add.expr = abline(h = c(0.5,5.5),v= c(0.5,5.5), lwd = 3),
# 
#           # cellnote=round(mat_data,digits=1),
#           notecex=1.5
# )
# mtext("Proportion vaccinated",side=1,cex=1.5,line=2.5,at=par("usr")[1]+.6*diff(par("usr")[1:2]))
# mtext("Proportion quarantined",side=4,cex=1.5,line=1,at=par("usr")[1]+0.4*diff(par("usr")[1:2]))
my.labs=c(0,10,20,30,40,50,60,70,80,90,100)
my.at=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
my.brks=seq(0, 100, by=10)
myColorkey <- list(at=my.brks, labels=list(at=my.brks, labels=my.labs,cex=1.5),height=1,width=3, space="top")
font.settings <- list(
  fontsize = 20,
  cex=1,
  fontfamily = "serif")
mapTheme <- c(rasterTheme(region=c("white","black"),par.xlab.text=font.settings,par.ylab.text=font.settings,axis.text=font.settings))  

x.scale <- list( alternating=1,cex=1.5)
y.scale <- list( alternating=1,cex=1.5)

levelplot(mat_data,
          xlab=list(label="Proportion isolated"),
          ylab=list(label="Proportion vaccinated"),
          par.settings=mapTheme,
          at=my.at,
         scales=list(x=x.scale, y=y.scale),
          colorkey=myColorkey,
          margin=F)
           

dev.off()


