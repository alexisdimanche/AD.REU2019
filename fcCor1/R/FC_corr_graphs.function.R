
####################################################################################
## plot.pheno.func
###################################################################################

#' Plot histogram of phenotype distribution for Both variables (Male and Female)
#'
#' Gives a phenotype vs.frequency histogram and overlays both variable. Ability to campare
#' distribution of the phenotype for both variables
#' @param pheno.mean.M data frame mean Male phenotype values for all phenotypes in dataset
#' @param pheno.mean.F data frame mean Female phenotype values for all phenotypes in dataset
#' @param pheno.of.interest specific name of phenotype of interest
#' @return p, a single histogram of the overlapping phenotype value distribution for both variables
#' @example
#' plot.pheno.func(pheno.mean.M = pheno.mean.M, pheno.mean.F=pheno.mean.F, pheno.of.interest="NMR_BF%_8wks", xlab = "BF % at 8wks")
#'
plot.pheno.func=function(pheno.mean.M=NULL, pheno.mean.F=NULL, pheno.of.interest="", xlab = ""){

  library(ggplot2)
  #Find where the pheno values are that correspond to the pheno of interest
  ind.pheno.m = which(rownames(pheno.mean.M) == pheno.of.interest)
  ind.pheno.f=which(rownames(pheno.mean.F)==pheno.of.interest)

  #extracts  mean pheno values for M and F for a specific "pheno.of.interest"
  pheno.interesting.M=pheno.mean.M[ind.pheno.m,]
  pheno.interesting.F=pheno.mean.F[ind.pheno.f,]

  #transpose data frame
  pheno.interesting.M.t=t(pheno.interesting.M) %>% as.data.frame()
  pheno.interesting.F.t=t(pheno.interesting.F) %>% as.data.frame()

  #Assign names (Female and Male) to values
  pheno.interesting.M.t$sex <- 'M'
  pheno.interesting.F.t$sex <- 'F'
  pheno.interesting.combined <- rbind(pheno.interesting.M.t, pheno.interesting.F.t)
  x=pheno.interesting.combined[1] %>% as.matrix()

  p <- ggplot(pheno.interesting.combined, aes(x, fill=sex))+geom_histogram()+scale_fill_manual("sex", values = c("F"="red", "M"="blue"))+xlab(xlab)+theme(text=element_text(size=18))
  return(p)
}



######################################################################################################
####### plot.hist.gene.func
######################################################################################

#' Plot  FC geno vs FC pheno scatter plot and histogram for a gene
#'
#' Function returns two plots 1)gene expression FC vs. phenotype FC for a specific gene and pheno
#' 2 )histogram FC for a specific gene
#' @param cors1 a data frame with the correlation values for a specific phenotype
#' @param dPhno a data frame with FC values for a specific phenotype for each strain
#' @param gene.probe the probe name associated to the gene of interest. Use maxcor.gene.name.func
#' @param pheno.of.interest specific name of phenotype of interest
#' @return two plots hist and plot
#' plot = scatter plot of gene expression FC vs. phenotype FC
#' hist = gene expression FC histogram
#' The correlation value of the gene of interest
#' @example
#' plot.hist.gene.func(cors1=cors1.NMR, dPhno=dPhno.NMR, gene.probe="1422582_at", pheno.of.interest="NMR_BF%_8wks")
#'
plot.hist.gene.func =function(cors1=NULL, dPhno=NULL, gene.probe="", pheno.of.interest="", ylab = "" ){

  #sort indices for specific gene
  ind.e = which(rownames(dExpr) == gene.probe)
  ind.p = which(rownames(dPhno) == pheno.of.interest)

  #Determine gene name for graph title
  gene = expr.formatted.gene.annotation$gene[which(expr.formatted.gene.annotation$probe == gene.probe)]

  #Compute Scatter plot and histogram
  xlab = paste(c(gene, "Expression FC (log2)"), collapse=" ")
  y = dPhno[ind.p,] %>% data.matrix %>% as.numeric
  x = dExpr[ind.e,] %>% data.matrix %>% as.numeric
  par(mar=c(5,5,2,2))
  plot = plot(x, y,ylab=ylab, xlab=xlab, cex.axis = 1.5, cex.lab = 2)

  abline(lm(y~x))

  print(plot)
  #hist = hist(dExpr[ind.e,] %>% data.matrix %>% as.numeric, xlab=gene, font = 2, font.lab = 8, font.axis = 4)
  #print(hist)

  r = bicorAndPvalue(dExpr[ind.e,] %>% data.matrix %>% as.numeric, dPhno[ind.p,] %>% data.matrix %>% as.numeric, use="pairwise.complete.obs")
  return(r$bicor)
}

#############################################################################################

###########################################################################################################
#####plot.hist.phenoFC.func
#############################################################################################

#' Histogram of the strain specific FC of a specific phenotype
#'
#' The function returns a histogram showing the distribution of the strain specific Fold-Change of a specific phenotype
#' @param dPhno.interest data frame with all the phenotype FC for each strain of a specific phenotype
#' @param dPhno.new.name Name that you want to assign to the reorganized data frame for the pheno of interest
#' @param xlab x axis name for histogram
#' @return one histogram and a list called 'resu' with the reorganized data frame 'DPhno.interest'
#' histogram of the Fc of a desired phenotype
#' reorganized data frame for a specific phenotype containing pheno value for each strain (NA's excluded)
#' @example
#' plot.hist.phenoFC.results.NMR = plot.hist.phenoFC.func(dPhno.interest= dPhno.NMR, xlab="NMR Body Fat % at 8 Weeks FC")
#'
plot.hist.phenoFC.func = function(dPhno.interest=NULL, xlab=""){

  str.ids = colnames(dPhno.interest)

  #switch orientation of data frame
  dPhno.interest = t(dPhno.interest)
  dPhno.interest = cbind(str.ids,dPhno.interest[,1]) %>% as.data.frame(stringsAsFactors=F)

  #rename colums and take out NA values
  names(dPhno.interest) = c("strain", "phenotype")
  dPhno.interest = na.omit(dPhno.interest)

  dPhno.interest$phenotype<-as.numeric(as.character(dPhno.interest$phenotype)) #for ggplot, geom_histogram to work

  #Plot histogram for phenotype fold change across strains
  p <- ggplot(dPhno.interest, aes(phenotype))+geom_histogram()+xlab(xlab)+theme(text=element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
  print(p)

  resu=list(dPhno.interest=dPhno.interest)
  return(resu)
}


##########################################################################################################
#########strain.spec.pheno.dec.func
#####################################################################

#' Barplot phenotype FC for every strain for a specific phenotype
#'
#' The function plots a color annoted barplot representing the FC for every strain for a specific phenotype.
#' Red represents a positive FC meaning the Female phenotype value is greater, Blue is the opposite
#' @param dPhno.interest specific phenotype FC values for every acceptable (no NA's) strains. Use the
#' output of the plot.hist.phenoFC.func)
#' @param gene.probe probe number of the gene of interest. Use gene with nax correlation. Refer to maxcor.gene.name.func
#' @param plot.title character of desired plot title
#' @return b , a barplot with the FC for every strain for a specific phenotype
#' @example
#' strain.spec.pheno.dec.func(dPhno.interest= dPhno.NMR.mod, gene.probe = "1422582_at", plot.title = "NMR BF % at 8wks")
#'
strain.spec.pheno.dec.func=function(dPhno.interest=NULL,gene.probe=NULL, plot.title =""){

  #organize specific gene and phenotype the same
  ind = which(rownames(dExpr) == gene.probe) # get expression FC from dExpr for a single gene
  df.gene = dExpr[ind,] %>% data.matrix #data frame with only the dExpr of a specific gene
  str.ids = colnames(df.gene) #strain ID's for every dExpr value
  df.gene = cbind(str.ids,df.gene[1,]) %>% as.data.frame(stringsAsFactors=F) #Bind strain and express value
  names(df.gene) = c("strain","expr") #assign names to columns

  #Organize strain names in the same order for gene and phenotype !Specific to each phenotype!
  pheno.strains = dPhno.interest$strain #After function having run plot.hist.phenoFC.func
  ind.g = sapply(pheno.strains,function(x){which(df.gene$strain==x)})
  df.gene = df.gene[ind.g,]
  all(dPhno.interest$strain == df.gene$strain)

  # sort indices for the phenotype. FC in decreasing order
  dPhno.interest$phenotype = dPhno.interest$phenotype %>% as.numeric
  dec.ind = order(dPhno.interest$phenotype,decreasing=T)

  #bar graph phenotype level for every strain
  pltdat.ph = dPhno.interest$phenotype[dec.ind] %>% as.numeric
  indpos = which(pltdat.ph > 0)
  indneg = which(pltdat.ph < 0)
  col = rep("red",length(c(indpos,indneg)))
  col[indneg] = "blue"
  b<-barplot(pltdat.ph, col=col, xlab="Strains", ylab="Fold-Change", main=plot.title, cex.axis = 1.5, cex.lab = 2)
  print(b)
}
######################################
##########################################################################################################
########strain.spec.gene.func
#####################################################################

#' Function to compute a strain specific gene expression Fold-Change barplot
#'
#' The function plots a color annoted barplot representing the FC for every strain for a specific gene.
#' Red represents a positive FC meaning the Female gene expression value is greater, Blue is the opposite
#' The order of the strain is the same as for the strain.spec.pheno.dec.func for the same dPhno.interest and gene
#' @param dPhno.interest specific phenotype FC values for every acceptable (no NA's) strains. Use the
#' output of the plot.hist.phenoFC.func
#' @param gene.probe probe number of the gene of interest. Use gene with nax correlation. Refer to maxcor.gene.name.func
#' @param plot.title character of desired plot title
#' @return b , a barplot with the FC for every strain for a specific phenotype
#' @example
#' strain.spec.gene.func(dPhno.interest = dPhno.NMR.mod, gene.probe = "1422582_at", plot.title = "Lep")
#'
strain.spec.gene.func = function(dPhno.interest = NULL, gene.probe = "", plot.title = ""){

  #organize specific gene and phenotype the same
  ind = which(rownames(dExpr) == gene.probe) # get expression FC from dExpr for a single gene
  df.gene = dExpr[ind,] %>% data.matrix #data frame with only the dExpr of a specific gene
  str.ids = colnames(df.gene) #strain ID's for every dExpr value
  df.gene = cbind(str.ids,df.gene[1,]) %>% as.data.frame(stringsAsFactors=F) #Bind strain and express value
  names(df.gene) = c("strain","expr") #assign names to columns

  #Organize strain names in the same order for gene and phenotype
  pheno.strains = dPhno.interest$strain #After function having run plot.hist.phenoFC.func
  ind.g = sapply(pheno.strains,function(x){which(df.gene$strain==x)})
  df.gene = df.gene[ind.g,]
  all(dPhno.interest$strain == df.gene$strain)

  # sort indices for phenotypes in decreasing order
  dPhno.interest$phenotype = dPhno.interest$phenotype %>% as.numeric
  dec.ind = order(dPhno.interest$phenotype,decreasing=T)

  #GENE level bar graph
  pltdat.gene = df.gene$expr[dec.ind] %>% as.numeric
  indpos = which(pltdat.gene > 0)
  indneg = which(pltdat.gene < 0)
  col = rep("red",length(c(indpos,indneg)))
  col[indneg] = "blue"
  b<-barplot(pltdat.gene, col=col, xlab="Strains",ylab="Expression FC (log2)", main = plot.title, font.main = 2, cex.axis=1.5, cex.lab=2)
  print(b)
}
######################################
##########################################################################################################
##plot.expr.func
###################################################################################

#' Plot histogram of the distribution of a specific gene for Both variables (Male and Female)
#'
#' Gives a gene expression vs.frequency histogram and overlays both variable. Ability to campare
#' distribution of the gene expression for both variables
#' @param expr.mean.M.unique data frame mean Male gene expression values for all genes in dataset
#' @param expr.mean.F.unique data frame mean Female gene expression values for all genes in dataset
#' @param gene.probe specific probe number for the gene of interest
#' @param gene.name name of the gene for the title
#' @return p, a single color coded histogram of the overlapping gene value distribution for both variables
#' @example
#' plot.expr.func(expr.mean.M.unique = expr.mean.M.unique, expr.mean.F.unique = expr.mean.F.unique, gene.probe = "1422582_at", gene.name = "Lep") #for NMR
#'
plot.expr.func=function(expr.mean.M.unique=NULL, expr.mean.F.unique=NULL, gene.probe = "", gene.name = ""){

  library(ggplot2)
  #extracts  mean expression values for M and F for a specific "gene.probe"
  ind.g.m = which(rownames(expr.mean.M.unique) == gene.probe)
  ind.g.f=which(rownames(expr.mean.F.unique) == gene.probe)
  gene.interesting.M=expr.mean.M.unique[ind.g.m,]
  gene.interesting.F=expr.mean.F.unique[ind.g.f,]
  gene.interesting.M.t=t(gene.interesting.M) %>% as.data.frame()
  gene.interesting.F.t=t(gene.interesting.F) %>% as.data.frame()
  gene.interesting.M.t$sex <- 'M'
  gene.interesting.F.t$sex <- 'F'
  gene.interesting.combined <- rbind(gene.interesting.M.t, gene.interesting.F.t)
  x=gene.interesting.combined[1] %>% as.matrix()

  p <- ggplot(gene.interesting.combined, aes(x, fill=sex))+geom_histogram()+scale_fill_manual("sex", values = c("F"="red", "M"="blue"))+xlab("Expression (log2)")+ylab("Frequency")+ggtitle(gene.name)+theme(text=element_text(size=18))
  return(p)
}

###########################
########################################################################################################
#####plot.hist.ExprFC.func
#############################################################################################

#' Histogram of the strain specific FC of a specific gene
#'
#' The function returns a histogram showing the distribution of the strain specific Fold-Change for a specific gene
#' @param gene.probe specific probe number for the gene of interest
#' @param gene.name name of the gene for the title
#' @return one histogram
#' histogram showing the distribution of the FC of a desired gene
#' @example
#' plot.hist.ExprFC.func(gene.probe = "1422582_at", gene.name = "Lep")
#'
plot.hist.ExprFC.func = function(gene.probe = "", gene.name = ""){


  ind.dExpr = which(rownames(dExpr) == gene.probe)
  dExpr.interest = dExpr[ind.dExpr,]
  dExpr.interest.t = t(dExpr.interest) %>% as.data.frame(stringsAsFactors=F)

  names(dExpr.interest.t) = c("expressionFC")
  dExpr.interest.t$expressionFC<-as.numeric(as.character(dExpr.interest.t$expressionFC))

  p <- ggplot(dExpr.interest.t, aes(expressionFC))+geom_histogram()+xlab("Expression Fold-Change (log2)")+ylab("Frequency")+ggtitle(gene.name)+theme(text=element_text(size=18), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16))
  print(p)
}


###############################################################
############ deg.ma.plot
##############################################################

#' function for plotting DEG analysis results
#'
#' The function outputs a pdF of a scatter plot of the average expression vs the log2 Fold_change for every gene
#' Only the genes satisfying threshold values for log2FC and pvalue are color annoted: red = positive FC, blue = negative FC
#' @param data is the "output" file from the function deg.analysis containing the FC for every gene and the AVG expression
#' across MAle and female (and other information of statistical significance)
#' @param fname file name of the pdf
#' @param fc.cut minimum absolute value for the FC threshold
#' @param fdr.cut maximum p value for color annotation
#' @return PDF containing a color annotated scatter plot of the Avg expression vs log2FC for every gene
#' @example
#' deg.ma.plot(data=output,fname="decode_maplot2.pdf")

deg.ma.plot = function(data=NULL,fname=NULL,fc.cut=2,fdr.cut=0.05){
  #format of PDF
  pdf(fname,height=4,width=4)

  #Select gene satisfying threshold conditions
  outf = output %>% filter(logFC > log2(fc.cut), adj.P.Val < fdr.cut)
  outm = output %>% filter(logFC < -log2(fc.cut), adj.P.Val < fdr.cut)

  #Plot scatter plot and assign colors
  plot(output$AveExpr, output$logFC, col="gray",
       xlab="Average expression",ylab="Log2 fold change")
  points(outf$AveExpr, outf$logFC, col="red")
  points(outm$AveExpr, outm$logFC, col="blue")
  abline(h=0,lty=2)
  dev.off()
}

#######################################################
######## comparecors.2methods.func
#########################################################

#' Comparing the distribution of correaltions for DEG's and FCcorG's for a specific phenotype
#'
#' @param DEG.cordata data frame with the correaltion of the genes that are most highly differentially
#' expressed. Use to select.def.func to compute
#' @param FCcorG data frame of correlation for most significant Fold Change genes use select.FCcorG.func
#' @return a boxplot graph with two boxplots comparing the correlation distribution for DEG's and FCcorG for a
#' specific phenotype
#' @example
#' comparecors.2methods.func(DEG.cordata = cors1.DEG.NMR, FCcorG.cordata = results.NMR.select.FCcorG.func, ylab="Correlation NMR BF%")
#'
comparecors.2methods.func = function(DEG.cordata = NULL, FCcorG.cordata = NULL, ylab = ""){

  FCcorG = data.frame(method = "FCcorG", correlation = abs(FCcorG.cordata$bicor))
  DEG = data.frame(method = "DEG", correlation = abs(DEG.cordata$bicor))

  boxplot.data = rbind(DEG, FCcorG)


  bp <- ggplot(boxplot.data, aes(x=method, y=correlation))+geom_boxplot()

  print(bp + My_Theme + ggtitle("Correlation for FCcorGenes and DEG's") +
          xlab("Methods") + ylab(ylab))
}

##############################################################
#####################venn.diagram.func
#######################################################

#' Function to create a Venn diagram between genes included
#' in the DEG analysis and FCcorG analysis for a specific phenotype
#'
#' @param DEG.cordata data frame with the correaltion of the genes that are most highly differentially
#' expressed. Use to select.def.func to compute
#' @param FCcorG.cordata data frame of correlation for most significant Fold Change genes use select.FCcorG.func
#' @param phenotype.of.interest name of the phenotype that we are looking at (used for the title of the diagram)
#' @example
#' venn.NMR.results = venn.diagram.func(DEG.cordata = cors1.DEG.NMR , FCcorG.cordata = results.NMR.select.FCcorG.func, phenotype.of.interest = "NMR phenotype")
#'
venn.diagram.func = function(DEG.cordata = NULL, FCcorG.cordata = NULL, phenotype.of.interest = ""){

  library(gridExtra)
  library(VennDiagram)

  #Detremine how many genes belong to each area of the Venn diagram
  area1 = nrow(DEG.cordata)
  area2 = nrow(FCcorG.cordata)

  #Determine the genes in the cross area (present in both DEG and FCcorG analysis)
  cross.area = DEG.cordata[DEG.cordata$probe %in% (FCcorG.cordata$probe),]
  cross.area.number = nrow(cross.area)

  g = draw.pairwise.venn(area1, area2, cross.area.number, euler.d = TRUE, category = c("DEG", "FCcorG"), fill = c("yellow", "purple"), col = c("yellow", "purple"), scaled = TRUE,
                         cat.pos = c(0,0) , cat.dist = rep(0.025, 2), cex = c(4,4,4), cat.cex = c(2, 2))
  #print(g)
  grid.arrange(gTree(children = g), top = phenotype.of.interest, bottom = "")

  return(cross.area)
  dev.off()

}
##############################################################################
##########heatmap.func
###########################################################

#' Creates a heatmap of the FC across strains
#'
#' Function that uses the result of select.FCcorG.func to create a heatmap representing the fold change across strains
#' for all FC genes for a specific phenotype
#' @param FCcorG.cordata data frame of correlation for most significant Fold Change genes, use select.FCcorG.func
#' @param dExpr dExpr data
#' @param PDF.name
#' @param title desired title of the heatmap
#' @param cellwidth desired width of cell. Has to be adjusted depending on number of strains and genes
#' @param cellheigth desired height of cell. Has to be adjusted depending on number of strains and genes
#' @return a PDF file with the heatmap
#' @return ah , a list  with the indices and other info about the heatmpa
#' @example
#' heatmap.NMR.results = heatmap.func(FCcorG.cordata = results.NMR.select.FCcorG.func, dExpr = dExpr, PDF.name = "heatmap.NMR", title = "NMR FCcorGenes", fc.cut.heat = 1, cellwidth = 0.2, cellheight = 3)
#'
heatmap.func = function(FCcorG.cordata = NULL, dExpr = dExpr, PDF.name = "", title = "", fc.cut.heat = NULL, cellwidth = 2, cellheight = 5){

  library(NMF)
  #Creates data frame of Fold Changes for only the FCcorGenes of interest
  ind.h = dExpr[(rownames(dExpr) %in% FCcorG.cordata$probe),]
  ind.h.t = t(ind.h) #transpose
  ind.h.t %>% data.matrix

  #Converts all values greater than to 1 to -1 and values smaller than -1 to -1
  #For visual clarity in the heatmap
  ind.h.t[ind.h.t < -fc.cut.heat] <- -fc.cut.heat
  ind.h.t[ind.h.t >fc.cut.heat] <- fc.cut.heat

  #Creates heatmap as pdf with personlizable cell size, axis name and title
  pdf(PDF.name, width=10, height=10)
  ah = aheatmap(ind.h.t ,breaks=0, Rowv = TRUE, Colv = TRUE, labRow = "Strains", labCol = "Genes", cellwidth = cellwidth, cellheight = cellheight
  )
  print(ah)
  dev.off()

  #returns list containing indices
  return(ah)
}

##########################################################################
#######pheno.heatmap. func
#################################################

#' Function to make a single column heatmap of the phenotype following specific indices
#'
#' @param results.heatmap the file corresponding to the output of the heatmap.func. It is specific to each phenotype and contains the indices of the strains
#' @param dPhno.pheno data frame of the phenotype FC for a specific phenotype across all strains
#' @param PDF.name name of the output PDF
#' @param cellwidth desired width of cell
#' @param cellheigth desired height of cell
#' @return a PDG with the single column heatmap and a file containing indices
#' @example
#' res.pheno.heatmap.NMR = pheno.heatmap.func(results.heatmap = heatmap.NMR.results, dPhno.pheno = dPhno.NMR, PDF.name = "heatmap.pheno.NMR" , cellwidth = 10, cellheight = 5)
#'
pheno.heatmap.func = function(results.heatmap = NULL, dPhno.pheno = NULL, PDF.name = "", cellwidth = 10, cellheight = 5 ){

  library(stats)
  library(data.table)

  ind.pheno = results.heatmap$rowInd #get indices from gene heatmap
  dPhno.heat = t(dPhno.pheno)
  dPhno.heat = dPhno.heat[ind.pheno,]  %>% as.data.frame(stringsAsFactors=F)
  dPhno.heat[is.na(dPhno.heat)] <- 0 #replace all NA by 0 values
  dPhno.heat = dPhno.heat %>% mutate(id = row_number()) #add column with indices
  names(dPhno.heat) = c("dPhno.FC", "order") #name columns


  pdf(PDF.name, width=10, height=10)
  ph = aheatmap(data.frame(dPhno.heat$dPhno.FC), Colv=NA,Rowv=NA, cellwidth = cellwidth, cellheight = cellheight)

  print(ph)
  dev.off()

  return(ph)
}

##############
