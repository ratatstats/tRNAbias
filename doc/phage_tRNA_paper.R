## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, message=FALSE, warning=FALSE)
library("tRNAbias")
library("scales")
library("png")
library("grid")

figlabel <- function(text, displ=0, cex=NULL, ...) {
 
  ds <- dev.size("in")
  # xy coordinates of device corners in user coordinates
  x <- grconvertX(c(0, ds[1]), from="in", to="user")
  y <- grconvertY(c(0, ds[2]), from="in", to="user")
 
  # account for the fragment of the device that 
  # the figure is using
  fig <- par("fig")
  dx <- (x[2] - x[1])
  dy <- (y[2] - y[1])
  x <- x[1] + dx * fig[1:2]
  y <- y[1] + dy * fig[3:4]

  sw <- 2*strwidth(text, cex=cex) * 60/100
  sh <- 2*strheight(text, cex=cex) * 60/100
 
  x1 <- x[1] + sw 
  y1 <- y[2] - sh
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
 
  text(x1, y1+displ, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

## ----genomeSummary, fig.height=2.7, fig.width=7, fig.cap="\\label{fig:genomeSummary}Phage 2.275.O carries 25 tRNA genes and is a large phage in both capsid size (120 nm) and genome size (348,911 bp). (A) Estimation of transcriptional units and their timing of expression. Top plot shows genome position of KEGG annotated genes, and bottom plot shows time to reach half the maximum expression of that gene (bottom). Blue and red bars indicate genes on the positive and negative strand, respectively. Early genes tend to be polymerases and sigma factors, while late genes tend to be structural proteins. More detailed annotations can be found in Supplementary Figure 1. (B) Electron microscopy image of phage 2.275.O."----

data(phage.cdstab)
data(keggAnnots)
data(rnaseqCoverageRatio)

matchKeggs = keggAnnots[phage.cdstab$geneName]

level3 = sapply(strsplit(unlist(sapply(matchKeggs[sapply(matchKeggs,length)>0],unique)),"\\."),function(x){x[3]})
levelsToPlot = rev(c(
	"09122 Translation",
	"09121 Transcription",
	"09108 Metabolism of cofactors and vitamins",
	"09104 Nucleotide metabolism",
	"09123 Folding, sorting and degradation",
	"09143 Cell growth and death",
	"09124 Replication and repair"
))

geneMidpoints = (phage.cdstab$pos1+phage.cdstab$pos2)/2
levelmidpoints = sapply(1:length(levelsToPlot),function(i){
		  	pointsToPlot = sapply(matchKeggs,function(x){any(grepl(levelsToPlot[i],x))})
			mean(geneMidpoints[pointsToPlot])	
		  })

phagePart = rnaseqCoverageRatio[grepl("^NVP",rownames(rnaseqCoverageRatio)),]
phagecenters = apply(phagePart[phage.cdstab$geneName,1:6],1,function(x){sum((1:length(x))*x)/sum(x)})

phageinit = apply(phagePart[phage.cdstab$geneName,1:6],1,function(x){
	       temp = gppred((1:length(x)),x/max(x),l=0.55,tau=0.5)
	       temp$x[min(which( (temp$x > 1) & (temp$mu > 0.5) ))]
     })

phageGeneExpressionMag = apply(phagePart[phage.cdstab$geneName,1:6],1,max)

toWrite = names(matchKeggs)
annotName = sapply(matchKeggs[toWrite],function(x){
			if(is.null(x)){return("")}
			y = strsplit(x,"\\.K?[0-9]+[[:blank:]]+")
			y = sapply(y,function(z){z[length(z)]})
			names(sort(table(y),decreasing=TRUE))[1]
		})

regionsToHighlight = list(c(52000,99800),c(217300,318000))

par(mar=c(0,0,0,0))
plot(0,0,col=0,ylim=c(-13,length(levelsToPlot)+4),
     xlim=padrange(c(phage.cdstab$pos1,phage.cdstab$pos2),c(0.39,0.42)),
     xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
geneOffset = 1
geneYPositions = (phageinit - max(phageinit))*3 - geneOffset
colorsToUse = rep(blues9[8],length(phageGeneExpressionMag))
colorsToUse[phage.cdstab$strand==-1] = reds9[8]

for(i in 1:length(regionsToHighlight)){
	eachRegion = regionsToHighlight[[i]]
	polygon(rep(eachRegion,each=2),c(padrange(geneYPositions),rev(padrange(geneYPositions))),col="#F0F0F0",border=NA)
	if(i==1){text(mean(eachRegion),min(geneYPositions)-1.5,"DNA Polymerase,\nRNA Polymerase Sigma Factors\nZinc Finger Protein",adj=c(0.5,1),cex=0.55)
	}else{text(mean(eachRegion),min(geneYPositions)-1.5,"Major Capsid Protein,\nTail Fiber, Baseplate,\nScaffolding Protein",adj=c(0.5,1),cex=0.55)}
}

segments(phage.cdstab$pos1,geneYPositions,phage.cdstab$pos2,geneYPositions,lwd=5,lend="butt",col=colorsToUse)
for(i in 1:length(levelsToPlot)){
	pointsToPlot = sapply(matchKeggs,function(x){any(grepl(levelsToPlot[i],x))})
	points(geneMidpoints[pointsToPlot],rep(i,sum(pointsToPlot)),pch=16,cex=0.7,col=maxDiffColors[i])
	segments(min(geneMidpoints),i,max(geneMidpoints),i,col=maxDiffColors[i],lwd=0.5)	
}
text(rep(min(c(phage.cdstab$pos1,phage.cdstab$pos2))-0.01*diff(range(c(phage.cdstab$pos1,phage.cdstab$pos2))),length(levelsToPlot)),
     1:length(levelsToPlot),gsub("^[0-9]+[[:blank:]]+","",levelsToPlot),cex=0.55,adj=c(1,0.5))
text(min(c(phage.cdstab$pos1,phage.cdstab$pos2))-0.23*diff(range(c(phage.cdstab$pos1,phage.cdstab$pos2))),mean(geneYPositions-0.5),"Transcription Timing",cex=0.7,srt=90)
toWrite = rep(TRUE,length(matchKeggs))
axis(2,at=((2:5) - max(phagecenters))*3 - geneOffset,
     labels=c(paste(seq(0,45,by=15)," min",sep="")),
     pos=(min(c(phage.cdstab$pos1,phage.cdstab$pos2))-0.01*diff(range(c(phage.cdstab$pos1,phage.cdstab$pos2)))),
     las=2,cex.axis=0.55,lwd=0.5)

trnaYPositions = rep(length(levelsToPlot)+1,dim(phage.trnatab)[1])
points((phage.trnatab[,"pos1"]+phage.trnatab[,"pos1"])/2,trnaYPositions,pch="|",cex=0.5,col=c("#000000",golds9[5])[as.numeric(phage.trnatab$intronstart>0)+1])
segments(min(geneMidpoints),trnaYPositions[1],max(geneMidpoints),trnaYPositions[1],lwd=0.5)	
text(min(c(phage.cdstab$pos1,phage.cdstab$pos2))-0.01*diff(range(c(phage.cdstab$pos1,phage.cdstab$pos2))),trnaYPositions[1],"tRNA",cex=0.55,adj=c(1,0.5))

axis(3, at = c(seq(1, 348911, by=100000),348911), pos = length(levelsToPlot)+2.5, cex.axis=0.55,lwd=0.5,tck=-0.02,labels=FALSE)
text(c(seq(0, 348911, by=100000),348911)+2000,rep(length(levelsToPlot)+3.8),
     paste(round(c(seq(0, 348911, by=100000),348911)/1000),"kb"),cex=0.55,adj=c(1,0.5))

phageEM = readPNG(system.file("extdata","em.png",package="tRNAbias"))
grid.raster(phageEM,x=unit(1, "npc"),y=unit(0.49, "npc"),height=unit(0.95, "npc"),just="right")

text(-150000,length(levelsToPlot)+2.8,"A",cex=1.2)
text(380000,length(levelsToPlot)+2.8,"B",cex=1.2,col="#FFFFFF")

## ----transcriptomeCodonBias, fig.height=3.2, fig.width=8, fig.cap="\\label{fig:transcriptomeCodonBias}Codon usage bias introduced by the phage tRNA pool is more pronounced in late genes than early genes. (A) A multidimensional scaling plot of phage and host proteins using Shannon-Jensen Divergence of the codon distributions shows that codon usage difference between phage and host. Points representing the codon recognition capacities of the tRNA pool for each organism are overlayed. Points representing the average codon usage for each organism are also overlayed. As the axes carry little meaning, they have been omitted. (B) Preference for the phage tRNA vs. the host tRNA pool (slant). (Zero signifies ambivalence.) (C) Mean codon usage for host and phage. (D) Slant toward the phage tRNA pool vs. timing of expression depicted as center of mass of RNA expression for the first round of infection. Note that this is different from the expression timing described in figure 1. The center of mass in this plot indicates how quickly the RNA transcript is degraded, while the time to half maximum expression shown in figure 1 summarizes transcription timing."----
data(smoothTranslate)
data(codonTable)
data(phageCodonBias)
data(hostCodonBias)
data(phageCodonProj)
data(hostCodonProj)
data(rnaseqCoverageRatio)
data(codonTable)
data(combinedTable.scale)
# data(phage.trnatab)

# Setting up the data:
trnaWithIntron = phage.trnatab$geneName[phage.trnatab$intronstart>0]

hostCodonTable = codonTable[grepl("^BCV",rownames(codonTable)),]
phageCodonTable = codonTable[grepl("^NVP",rownames(codonTable)),]

hostCodonCent = apply(hostCodonTable,2,mean)
phageCodonCent = apply(phageCodonTable,2,mean)

phagePart = rnaseqCoverageRatio[grep("^NVP|^phage",rownames(rnaseqCoverageRatio)),]
hostPart = rnaseqCoverageRatio[grep("^BCV|^host",rownames(rnaseqCoverageRatio)),]

phagecenters = apply(phagePart[,1:6],1,function(x){sum(((1:length(x)) )*x)/sum(x)})
phagecenters = phagecenters[names(phageCodonBias)]

phagetRNAPart = phagePart[grep("^phage",rownames(phagePart)),]
phagetRNAPart = phagetRNAPart[!(rownames(phagetRNAPart) %in% trnaWithIntron),]
phagetRNAcenters = apply(phagetRNAPart[,1:6],1,function(x){sum(((1:length(x)) )*x)/sum(x)})

phage.loess = loess(phageCodonBias ~ phagecenters)
phage.xvals = seq(min(phagecenters), max(phagecenters), length.out=200)
phage.loesspredict = predict(phage.loess, data.frame(phagecenters = phage.xvals), se = TRUE)

hosttRNAPart = hostPart[grep("^host",rownames(hostPart)),]
hosttRNAlows = apply(hosttRNAPart[,1:6],1,which.min) 


layout(matrix(c(4,4,1,2,3,3),nrow=2),widths=c(1.7,1.73,2))

x = c(rep(1,dim(hostCodonTable)[1]),rep(2,dim(phageCodonTable)[1]))
y1 = c(hostCodonProj,phageCodonProj)
y2 = c(hostCodonBias,phageCodonBias)

par(mar=c(4.2,3.7,2.8,1))
scatterBox(x,y2,ptcols=c(purples9[7],golds9[7])[x],
	   outline=FALSE,nbreaks=30,
	   border=c(purples9[9],golds9[9]),
	   names=c("",""), xlab="",
	   horizontal=TRUE,
	   col=alpha("#FFFFFF",0))
title(xlab="host <- tRNA Preference -> phage",line=2.5)
title(main="tRNA Pool Preference",font.main=1)
axis(2,at=c(1,2),labels=c("Host\nGenes","Phage\nGenes"),lty=0,las=2)
figlabel("B",cex=1.5)

scatterBox(x,y1,ptcols=c(purples9[7],golds9[7])[x],
	   outline=FALSE,nbreaks=30,
	   border=c(purples9[9],golds9[9]),
	   names=c("",""), xlab="",
	   font.main=1,horizontal=TRUE,
	   col=alpha("#FFFFFF",0))
title(xlab="host <- Codon Similarity -> phage",line=2.5)
title(main="Codon Usage Similarity",font.main=1)
axis(2,at=c(1,2),labels=c("Host\nGenes","Phage\nGenes"),lty=0,las=2)
figlabel("C",cex=1.5)

par(mar=c(4.2,4.7,2.8,1))
plot(phagecenters, phageCodonBias, col=0, xlab="", ylab="",xaxt="n",
     main="Phage Genes by Timing", font.main=1,
     ylim=range(c(phageCodonBias)))#,hostCodonBias
abline(v=jitter(hosttRNAlows),col=purples9[3],lty=3)
abline(v=phagetRNAcenters,col=golds9[3])
points(phagecenters, phageCodonBias, col=golds9[7])
title(xlab="early <- Expression Timing -> late",line=2.5)
title(ylab="host <- tRNA Pool Preference -> phage",line=2.5)
lines(phage.xvals,phage.loesspredict$fit,lty=1,col=golds9[9],lwd=2)
lines(phage.xvals,phage.loesspredict$fit+1.98*sqrt(phage.loesspredict$se.fit^2),lty=2,col=golds9[9],lwd=1.2) 
lines(phage.xvals,phage.loesspredict$fit-1.98*sqrt(phage.loesspredict$se.fit^2),lty=2,col=golds9[9],lwd=1.2)
axis(1,at=1:10,labels=c("Preinf",paste(seq(0,120,by=15)," min",sep="")))
figlabel("D",cex=1.5)

par(mar=c(1,1,2.8,0))
plot(combinedTable.scale[-c(1,2,3,4),],main="Codon Usage MDS",
     col=c(colorRampPalette(c("#FFFFFF",purples9[8]))(6)[3],
           colorRampPalette(c("#FFFFFF",golds9[8]))(6)[3])[
             c(rep(1,dim(hostCodonTable)[1]),rep(2,dim(phageCodonTable)[1]))],
     xlim=padrange(combinedTable.scale[,1],0.06),
     ylim=padrange(combinedTable.scale[,2],c(0.06,0)),
     xaxt="n",yaxt="n",bty="n",font.main=1)
legend("bottomright",pch=15,
       col=c(purples9[8],golds9[8],"#FFFFFF"),
	legend=c("Host        ","Phage        ",""),pt.cex=1.8,bty="n")
title(xlab="MDS Axis 1",line=-0.7)
title(ylab="MDS Axis 2",line=-0)
#title(main="Codon Usage MDS",line=4)
points(combinedTable.scale[1:4,],col=rep(c(purples9[8],golds9[8]),2),cex=1.5,pch=16)
text(combinedTable.scale[1:4,1]+0.01,combinedTable.scale[1:4,2],labels=rep(c("tRNA","codon"),each=2),
       col=rep(c(purples9[9],golds9[9]),2),cex=1,pch=16,adj=c(0,0.5))
figlabel("A",cex=1.5)

## ----genomicCodonBias, fig.height=4.8, fig.width=7.3, fig.cap="\\label{fig:genomicCodonBias}Analysis of the phage and host genome supports the codon usage hypothesis. (A) tRNA content in the genomes of the phage and host. Less saturated bars indicate putative tRNA with introns. (B) Differences between the codon usages of phage and host. Darker-hued bars indicate codons that cannot be recognized by phage tRNA, given the wobble rules summarized by **dos Reis, et al.**"----

data(codonKey)
data(host.cdsseqs)
data(phage.cdsseqs)
data(host.trnatab)
data(phage.trnatab)
data(codonTable)
data(translateTable)


hostTable = codonTable[grepl("^BCV",rownames(codonTable)),]
phageTable = codonTable[grepl("^NVP",rownames(codonTable)),]

tRNASummary = table(c(host.trnatab$organism,phage.trnatab$organism),
		    factor(gsub("^[^-]+-([^-]+)$","\\1",c(host.trnatab$annot,phage.trnatab$annot)),
			   levels=rc(codonKey[,1]),ordered=TRUE),
		    c(host.trnatab$intronstart,phage.trnatab$intronstart)>0)
codonSummary = rbind(apply(hostTable,2,sum),apply(phageTable,2,sum))

codonSummProp = sweep(codonSummary,1,apply(codonSummary,1,sum),"/")

orcodon = (codonSummProp[2,]/(1-codonSummProp[2,])) / (codonSummProp[1,]/(1-codonSummProp[1,]))

phageTranslationCapacity = dosReisWobble %*% table(factor(gsub("^.*-","",phage.trnatab$annot[!phage.trnatab$intronstart>0]),levels = colnames(dosReisWobble))) > 0
hostTranslationCapacity = dosReisWobble %*% table(factor(gsub("^.*-","",host.trnatab$annot[!host.trnatab$intronstart>0]),levels = colnames(dosReisWobble))) > 0
# phageTranslationCapacity = translateTable %*% table(factor(phage.trnatab$geneName,levels = colnames(translateTable))) > 0
# hostTranslationCapacity = translateTable %*% table(factor(host.trnatab$geneName,levels = colnames(translateTable))) > 0

trnabarplot = function(x,...){
bartrick = matrix(0,2*dim(x)[2],dim(x)[2])
bartrick[cbind(1:dim(bartrick)[1],rep(1:dim(x)[2],each=2))] = x
barplot(bartrick,
	col=as.vector(rbind(aacolors[6,codonKey[,3]],aacolors[3,codonKey[,3]])),
	ylab="Count",las=2,yaxt="n",xaxt="n",font.main=1, ...)
axis(side=2,at=c(0,max(apply(x,2,sum))),cex=1)
}

aacolors = cbind(golds9,seagreens9,reds9,blues9,rep("#D3D3D3",9))
colnames(aacolors) = c("nonpolar","polar","acidic","basic","stop")

cexaa = 0.7

layout(matrix(c(1,2,3),3,1),heights=c(0.43,0.55,1.5),widths=c(1.8))

par(mar=c(1,4.9,1.8,0))
trnabarplot(t(tRNASummary["phage",,]),main="tRNA, Phage",cex.axis=cexaa,ylim=c(0,2.3))
figlabel("A",displ=0.3,cex=1.5)

par(mar=c(1,4.9,1.3,0))
trnabarplot(t(tRNASummary["host",,]),main="tRNA, Host",cex.axis=cexaa)
axis(side=1,gsub("\\([^)]+\\)","",paste(codonKey[,2],gsub("T","U",rc(codonKey[,1])),sep="-")),
     at=1:dim(codonKey)[1]*1.2-0.5,cex.axis=cexaa,cex=cexaa,las=2,tick=FALSE,pos=1.5)

par(mar=c(0,4.9,2.5,0))
x = log2(orcodon)
barplot(x,col=diag(aacolors[6+((!phageTranslationCapacity) & hostTranslationCapacity)*3,codonKey[,3]]),
	ylim=padrange(x,c(0.27,0.25)),xaxt="n", ylab="Host  <-  log2( Odds Ratio )  ->  Phage")
text(1:length(x)*1.2-0.5,
	x + 0.12*diff(range(x,na.rm=TRUE))*((as.numeric(x>0)*2-1)),
	labels=paste(gsub("\\([^)]+\\)","",codonKey[1:length(x),2]),
		     gsub("T","U",codonKey[1:length(x),1]),sep="-"),srt=90,cex=cexaa)
title(main="Codon Odds Ratio, Phage to Host",line=-1.5, font.main=1)
legend("bottomright",pch=15,col=aacolors[6,],
	legend=colnames(aacolors),ncol=5,pt.cex=1.8,bty="n")
figlabel("B",displ=-0.55,cex=1.5)

## ----phageLifecycle, fig.height=6, fig.width=3.3, fig.cap="\\label{fig:phageLifecycle}Time course of phage 2.275.O. infection. (A) Phage plaques begin to appear 30 minutes into infection, and peaks at approximately 90 minutes. (B) Genomic qPCR results show that the host genome is degraded rapidly upon infection. Approximately a tenth of the hosts remain uninfected and begin to grow again. Then we can see a second round of infection starting at approximately 90 minutes. (C) The tRNA subset of RNA-seq shows that host tRNA (and generally, most host RNA, with the exception of stress response genes) are degraded upon infection while phage tRNA rapidly increase. Reads are normalized to a firefly luciferase spike-in for each sample, as opposed to the total read count per sample."----
data(standardCurvesData)
data(qpcrData)
data(pfuData)
data(rnaseqCoverageRatio)

tRNAPart = rnaseqCoverageRatio[grep("^host|^phage",rownames(rnaseqCoverageRatio)),]

standardCurves.lm = lm(conc ~ Cq*primer, data=standardCurvesData)

ratioAmt = rep(1,length(qpcrData$sampleID))
ratioAmt[qpcrData$sampleID==1] = 2*0.36
# The preinfection sample is 2.78x more concentrated than other samples,
# and at the time of the experiment was diluted 2 fold for qPCR

allPrimers = c("host GroEL","host CTP Sythetase","phage GroEL","phage Major Capsid")
colorsToUse = c("#218EBF", "#9242C6", "#C93A62", golds9[6])
names(colorsToUse) = allPrimers

# Calculate the concentrations of DNA using the standard curve fits
calcConcentrations = predict(standardCurves.lm,qpcrData)
calcConcentrations = calcConcentrations + log10(ratioAmt * 5)
# Multiply by 5 to get into the same units as plaquing data
# calcConcentrations[primer=="phage GroEL"] = NA

outterxlim = c(1,dim(pfuData)[1]+3)

### Plot filtered and unfiltered phage

layout(matrix(1:3,nrow=3),heights=c(1.1,1,1.1))
par(mar=c(0.1,6.3,4,1))

plot(0, 0, col = 0,
     xlim=outterxlim,
     ylim=range(unlist(pfuData[,1])),
     xlab="",ylab="log10 ( Plaques / 5uL )",
     xaxt="n",bty="n",
     main="Infection Timecourse Characteristics",font.main=1)
abline(v = 6, col = golds9[2], lwd = 5)
lines(1:dim(pfuData)[1] + 1, pfuData[,1],
      col=colorRampPalette(c(golds9[6],"#FFCC33"))(3)[2],
      lwd = 1.5)
# lines(1:dim(pfuData)[1] + 1,pfuData[,2],lty=2,col=golds9[9],lwd=2)
# legend("bottomright",lty=c(1,2),col=golds9[c(6,9)],
#        legend=c("unfiltered pfu", "filtered pfu"),bty="n",cex=0.9,pt.cex=1,lwd=2)
figlabel("A",disp=-0.5,cex=1.5)

par(mar=c(0.1,6.3,0.1,1))
plot(qpcrData$sampleID,
     calcConcentrations,
     type = "n",
     xlab="",ylab="log10 ( qPCR copy number / 5uL )",
     xaxt="n",bty="n",yaxt="n",
     ylim=padrange(calcConcentrations,c(0.1,0.07)),
     xlim = outterxlim)
# polygon(c(0,0,11,11),c(0,10,10,0),col="#F8F8F8",border=NA)
abline(v = 6, col = golds9[2], lwd = 5)
axis(2,at=5:8)
points(qpcrData$sampleID,
     calcConcentrations,
     col=colorsToUse[qpcrData$primer],pch=16)
legend("bottomright",pch=16,col=colorsToUse,
       legend=names(colorsToUse),bty="n",cex=0.9,pt.cex=1)
figlabel("B",cex=1.5)

par(mar=c(4,6.3,0.1,1))
plot(0,0,
     xlim=outterxlim,#c(1,dim(tRNAPart)[2]),
     ylim=padrange(tRNAPart,c(0,0.06)),
     bty="n",xaxt="n",
     xlab="",ylab="tRNA Expression (Relative)")
abline(v = 6, col = golds9[2], lwd = 5)
for(x in 1:dim(tRNAPart)[1]){
  lines((1:dim(tRNAPart)[2]),tRNAPart[x,],lwd=1.5,
	col=c(purples9[8],colorRampPalette(c(golds9[6],"#FFCC33"))(3)[2])[
	  grepl("^phage",rownames(tRNAPart)[x])+1])
  }
legend("topright",col=c(purples9[8],colorRampPalette(c(golds9[6],"#FFCC33"))(3)[2]),legend=c("host tRNA","phage tRNA"),
       lwd=1.5,lty=1,bty="n",cex=0.9)
axis(1,at=1:10,labels=c("Preinf",paste(seq(0,120,by=15)," min",sep="")),las=2,srt=45)
figlabel("C",cex=1.5)

## ----diversityAnalysis, fig.height=3, fig.width=5.5, fig.cap="\\label{fig:diversityAnalysis}tRNA carried by the phage may supplement the degrading pool of host tRNA. (A) The probability of selecting, uniformly at random from the host genome, a tRNA array that is able to encode as many anticodons as that carried by the phage is 0.0016. In addition, contiguous stretches of tRNA in the host genome, which are typically thought to be the result of duplication events, encode very lowly diverse anticodons. The phage’s tRNA collection, therefore, appears to be the result of multiple acquisition and selection events. (B) Of the tRNA not carried by the phage, most are highly expressed by the host, and others correspond to codons not very highly used by the phage genome."----

### tRNA Diversity ###

hasintron = phage.trnatab$intronstart>0
B = 10000

set.seed(1000)
# Sample codons without replacement
sampleCodonsWithReplacement = function(){
	sample(host.trnatab$annot,sum(!hasintron),replace=TRUE)
}
repSampleCodons = replicate(B,sampleCodonsWithReplacement())

randDiversityCodons = apply(repSampleCodons,2,function(x){entropy(table(x)/length(x))})
phageDiversityCodons = entropy(table(phage.trnatab$annot[!hasintron])/sum(!hasintron))

# Are contiguous tRNA more or less unique?
indexhack = tapply(1:dim(host.trnatab)[1],host.trnatab$contig,function(x)x)
hostadjtRNA = rep(0,dim(host.trnatab)[1])
temp = tapply(host.trnatab$pos1,host.trnatab$contig,function(x){
		midinds = 1:length(x)
		sortorder = order(x)
		hasbreak = c(0,cumsum(diff(x[sortorder])>25000))
		adjtRNA = rep(0,length(hasbreak))
		adjtRNA[order(x)] = hasbreak
		adjtRNA
	})
hostadjtRNA[unlist(indexhack)] = unlist(temp)

hostContigCodonDiversity = tapply(host.trnatab$annot,paste(host.trnatab$contig,hostadjtRNA,sep=" "),
				function(x){entropy(table(x)/length(x))})


### And other characteristics ###

# P(codon) = sum_g P(codon | gene = g) * P(gene = g)

codonProps = sweep(codonTable,1,apply(codonTable,1,sum),"/")

hostPCodonGGene = codonProps[grepl("^BCV",rownames(codonProps)),]
phagePCodonGGene = codonProps[grepl("^NVP",rownames(codonProps)),]

hostPGene = rnaseqCoverageRatio[rownames(hostPCodonGGene),1]
phagePGene = apply(rnaseqCoverageRatio[rownames(phagePCodonGGene),2:6],1,mean)
hostPGene = hostPGene/sum(hostPGene)
phagePGene = phagePGene/sum(phagePGene)

hostPCodon = (t(hostPGene) %*% hostPCodonGGene)
phagePCodon = (t(phagePGene) %*% phagePCodonGGene)

# theoretical tRNA expression that would match codon usage
# P(trna) = sum_c P(trna | codon = c) * P(codon = c)

hosttRNAbinding = bindingAffinities[,grep("^host",colnames(bindingAffinities))]
tRNAtheorBias = sweep(hosttRNAbinding,1,apply(hosttRNAbinding,1,sum),"/")

hosttRNABias = hostPCodon %*% tRNAtheorBias
phagetRNABias = phagePCodon %*% tRNAtheorBias

hosttRNAExpression = tRNASamples[grep("^host",rownames(tRNASamples)),c("C0","C15","C30","C45","C60")]
hosttRNAExpression = sweep(hosttRNAExpression, 2, apply(hosttRNAExpression, 2, sum), "/")
hosttRNAExpression = apply(hosttRNAExpression,1,mean)

tRNACarriedByPhage = (gsub("^.*-(.*-.*)-.*$","\\1",names(hosttRNAExpression)) %in% phage.trnatab$annot[phage.trnatab$intronstart==0] + 1)

### Plot Code ###

layout(matrix(1:3,nrow=1),widths=c(1,1,0.35))

x = c(rep(2,length(hostContigCodonDiversity)),rep(1,length(randDiversityCodons)))
y = c(hostContigCodonDiversity,randDiversityCodons)

par(mar=c(5,4.75,4,1.25))
scatterBox(x,y,ptcols=c("#888888",purples9[7])[x],
	   outline=FALSE,pch=16,nbreaks=30,
	   border=c("#0A0A0A",purples9[9]),
	   names=c("",""), ylab="Anticodon Entropy",
	   col=alpha("#FFFFFF",0))
title(main="Phage tRNA\nAnticodon Diversity", font.main=1)
axis(1,at=c(2,1),labels=c("Contiguous\nHost tRNA","Randomly\nAcquired"),line=0.9,lty=0)#,cex.axis=0.8
abline(h=phageDiversityCodons,col=golds9[9],lwd=1.5,lty=2)
text(2.5,phageDiversityCodons-0.15,
     paste("Phage tRNA\npval = ",sum(randDiversityCodons >= phageDiversityCodons)/B),
     adj=c(1,1),col=golds9[9])
figlabel("A",cex=1.5)

phagebreaks = cut(phagetRNABias[,names(hosttRNAExpression)],breaks=9)
phagecolors = colorRampPalette(golds9[-1])(length(levels(phagebreaks)))

par(mar=c(5,4.75,4,1.25))
scatterBox(c("Absent\nfrom Phage","Carried\nby Phage")[tRNACarriedByPhage], hosttRNAExpression , 
	   ptcol = phagecolors[phagebreaks],
	   outline=FALSE,pch=16,las=2,yaxt="n",
	   ylab="Host tRNA Expression",names=c("",""),
	   col=alpha("#FFFFFF",0))
title(main="Host tRNA\nCharacteristics",font.main=1)
axis(1,at=c(2,1),labels=c("Has Phage\nAnalog","No Phage\nAnalog"),line=0.9,lty=0)
figlabel("B",cex=1.5)

legendtext = rep("",length(levels(phagebreaks)))
legendtext[2:7] = c("  high",expression("   "%up%""),
		  "Phage Codon","Recognition",
		  expression("   "%down%""),"  low")
par(mar=c(4,0,2.8,0))
plot(0,0,col=0,bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
legend("left",col=rev(phagecolors),bty="n",
       legend = legendtext,#legend=levels(phagebreaks),
       pch=15,cex=0.9,pt.cex=2.5)

