
##########################################
###expr.format.func
###########################################

#'Function Sets expression matrix and annotation
#'
#'@param expr0 which is the output from the read.table function
#'@returm data frame witht he formatted exopression data
#' @example
#' expr.format.results=expr.format.func(expr0)
#'
expr.format.func=function(expr0=NULL){
  #set expression matrix and annotation
  #takes about 15 min
  expr1 = reshape(expr0[,c(1,2,5)], idvar="probesetID", timevar="mouse_number", direction="wide")
  ## applies function to column names that gets rid of expression_value. and leaves just the number
  expr.mouse.nums = sapply(names(expr1)[2:ncol(expr1)],function(x){strsplit(x,"expression_value.")[[1]][2]})
  ## get annotation for strain values
  expr.mouse.strains = var.annot(vecdat=expr.mouse.nums, annot=annotation, invar="mouse_number", outvar="Strain")
  expr.mouse.sex = var.annot(vecdat=expr.mouse.nums, annot=annotation, invar="mouse_number", outvar="Sex")
  ## just expression values
  expr2 = expr1[,2:ncol(expr1)]
  # adds row names that are the probesetIDs
  rownames(expr2) = expr1$probesetID
  #rename columns with just the mouse number
  names(expr2) = expr.mouse.nums
  print("have done names")
  #group together the annotations created with the function previously
  #there is one annotation that is mouse number, strain, and sex
  expr2.samp.annotation = cbind(expr.mouse.nums, expr.mouse.strains, expr.mouse.sex) %>% as.data.frame(stringsAsFactors=F)
  names(expr2.samp.annotation) = c("mouse_number","strain","sex")
  print("right before slow step")
  #there is another annotation that is probesetID with gene symbol
  expr.genes = var.annot(vecdat=rownames(expr2), annot=expr0, invar="probesetID", outvar="gene_symbol")
  print("finished slow step")
  expr2.gene.annotation = cbind(rownames(expr2), expr.genes) %>% as.data.frame(stringsAsFactors=F)
  names(expr2.gene.annotation) = c("probe","gene")
  print("right before out")
  out = list(expr=expr2, expr.samp.annot=expr2.samp.annotation,expr.gene.annot=expr2.gene.annotation)
  return(out)
}

##########################################
###pheno.format.func
###########################################

#' The function sets phenotype matrix and annotation
#' Similar function than expr.format.func
#'
#' @param pheno
#' @return formatted phenotype data
#' @examples
#' pheno.format.results=pheno.format.func(pheno0)
#'
pheno.format.func=function(pheno=NULL){
  pheno1=reshape(pheno[,c(1,2,5)], idvar="trait_name", timevar="mouse_number", direction="wide")
  pheno.mouse.nums = sapply(names(pheno1)[2:ncol(pheno1)],function(x){strsplit(x,"value.")[[1]][2]})
  pheno.mouse.strains = var.annot(vecdat=pheno.mouse.nums, annot=annotation, invar="mouse_number", outvar="Strain")
  pheno.mouse.sex = var.annot(vecdat=pheno.mouse.nums, annot=annotation, invar="mouse_number", outvar="Sex")
  # phenotype expressions
  pheno2 = pheno1[,2:ncol(pheno1)]
  rownames(pheno2) = pheno1$trait_name
  names(pheno2) = pheno.mouse.nums
  #create an annotation matrix with mouse number, strain, and sex
  pheno2.samp.annotation = cbind(pheno.mouse.nums, pheno.mouse.strains, pheno.mouse.sex) %>% as.data.frame(stringsAsFactors=F)
  names(pheno2.samp.annotation) = c("mouse_number","strain","sex")
  out = list(pheno=pheno2, pheno.samp.annotation=pheno2.samp.annotation)
  return(out)
}

#########################################################
##filtering phenotype data
########################################################

#' Function to filter NA values from phenotype data
#'
#' @param pheno pheno.formatted data frame, output from pheno.format.func$pheno
#' @param na.prop threshold value
#' @examples
#' pheno6=remove.na.func(pheno=pheno.formatted, na.prop=0.3)
remove.na.func= function(pheno=NULL, na.prop=0.3){
  pheno=t(pheno)
  res=c()
  column.names=c()
  for (ii in 1:ncol(pheno)){
    x=pheno[,ii] %>% data.matrix() %>% as.numeric()
    # if(ii==3){break}
    if(length(which(is.na(x)==TRUE)) < na.prop*nrow(pheno)){
      res=cbind(res, pheno[,ii])
      new.column.name=colnames(pheno)[ii]
      column.names=c(column.names, new.column.name)
      column.names=t(column.names)
    }
  }
  colnames(res)=column.names
  pheno6 <- res


  return(pheno6)
}

############################################################
########strains.for.analysis.func
############################################################

#' selects strains for analysis
#'
#' @param expr2
#' @param pheno2
#' @param expr2.samp.annotation
#' @param pheno2.samp.annotation
#' @examples
#' strains.for.analysis.results=strains.for.analysis.func(expr2=expr.formatted, pheno2=t(pheno6), expr2.samp.annotation=expr.formatted.samp.annotation, pheno2.samp.annotation=pheno.formatted.samp.annotation)
#'
strains.for.analysis.func= function(expr2=NULL,pheno2=NULL,expr2.samp.annotation=NULL,pheno2.samp.annotation=NULL){

  # identify strains with both M and F - expression and phenotype
  #for the phenotype and expression data, select the strains corresponding to male and females
  # remove duplicates
  #using intersection: find a strain that has a phenotype associated with it and has male and female
  strainsM.expr = expr2.samp.annotation %>% filter(sex == "M") %>% select(strain) %>% unique
  strainsM.pheno = pheno2.samp.annotation %>% filter(sex == "M") %>% select(strain) %>% unique
  strainsF.expr = expr2.samp.annotation %>% filter(sex == "F") %>% select(strain) %>% unique
  strainsF.pheno = pheno2.samp.annotation %>% filter(sex == "F") %>% select(strain) %>% unique
  strainsMF = Reduce(intersect, list(strainsM.expr$strain,strainsM.pheno$strain,strainsF.expr$strain,strainsF.pheno$strain))

  # filter expression data for intersect strains
  # main objects:
  ##### expr3: gene expression matrix
  ##### expr3.samp.annotation: column/sample annotation for expr3; mouse, strain, sex
  ##### expr2.gene.annotation: row/probe annotation for expr3; probe, gene name

  #returns true or false for each annotation saying whether or not the strain has male and females and phenotype data
  expr.strain.inds = expr2.samp.annotation$strain %in% strainsMF

  #expr 3 annotations only include strains with both male and female
  expr3.samp.annotation = expr2.samp.annotation[expr.strain.inds,]
  #exressions for strains with M and F and unique
  expr3 = expr2[,expr.strain.inds]

  # do same thing you did for expr. selecting the ones that are in MF intersect
  # main objects:
  ##### pheno3: phenotype matrix
  ##### pheno3.samp.annotation: column/sample annotation for pheno3; mouse, strain, sex
  pheno.strain.inds = pheno2.samp.annotation$strain %in% strainsMF
  pheno3.samp.annotation = pheno2.samp.annotation[pheno.strain.inds,]
  ##this is where error is
  pheno3 = pheno2[,pheno.strain.inds]

  # modify organization
  # expr and pheno data are organized identically by mouse number
  # main objects:
  ##### pheno4: phenotype matrix
  ##### pheno4.samp.annotation: column/sample annotation for pheno4; mouse, strain, sex
  expr.nums = expr3.samp.annotation$mouse_number
  #gives the indices where the mouse number of pheno is equal to the mouse number from expr
  pheno.inds = sapply(expr.nums,function(x){which(pheno3.samp.annotation$mouse_number==x)}) %>% unlist
  # edits pheno3 annotation to only include annotations for the mouse numbers we have expressions for
  pheno4.samp.annotation = pheno3.samp.annotation[pheno.inds,]
  #same thing as above but for phenotype values
  pheno4 = pheno3[,pheno.inds]
  out=list(expr.samp.annot=expr3.samp.annotation, expr=expr3, pheno.samp.annot=pheno4.samp.annotation, pheno=pheno4)
  return(out)
}
############################################################
##############separate.by.sex.func
############################################################

#' The Function seperate.by.sex.func uses finds indices corresponding to Male and Female data and then seperates annotation,
#' expression, and phenotype data by sex
#'
#' @param samp.annotation data frame containing the strain and sex for every mouse number
#' @param expr  data frame containing gene expression data for each mouse number for each probe number
#' @param pheno data frame containing the phenotypic data for all recorded phenotypes across all mouse numbers
#' @return a list containing annotation, expression and phenotype data seperated by sex
#' @example
#' separate.by.sex.results=separate.by.sex.func(samp.annotation, expr.for.analysis,pheno.for.analysis)
#'
separate.by.sex.func=function(samp.annotation=NULL, expr=NULL, pheno=NULL){

  #locate which data corresponds to Male and Female
  ind.m = which(samp.annotation$sex=="M")
  ind.f = which(samp.annotation$sex=="F")

  anno.M = samp.annotation[ind.m,] # annotation, male
  anno.F = samp.annotation[ind.f,] # annotation, female

  #seperate expression data by sex
  expr.M = expr[,ind.m] # expression, male
  expr.F = expr[,ind.f] # expression, female

  #seperate phenotype data by sex
  pheno.M = pheno[,ind.m] # phenotype, male
  pheno.F = pheno[,ind.f] # phenotype, female

  out=list(anno.M=anno.M, anno.F=anno.F, expr.M=expr.M, expr.F=expr.F, pheno.M=pheno.M, pheno.F=pheno.F)
  return(out)
}


###################################################################
##########get.means
###########################################

#' function to get means and revised annotations
#'
#' gets the sex specific mean expression and mean phenotype data for each strain
#' @param expr data (can be pheno or expression data) specific to one sex
#' @param anno sex specific annotation data containing the strain for each mouse number
#' @param strains strain names present for a specific sex
#' @return list composed of 2 data frames expr & anno
#' @example
#' expr.M.mean = get.means(expr=expr.M, anno=anno.M, strains=unique(anno.M$strain))
#'
get.means = function(expr=NULL,anno=NULL,strains=NULL){
  mat = matrix(0, nrow(expr), length(strains)) %>% as.data.frame
  nreps = rep(0, length(strains))
  for(ii in 1:length(strains)){
    ind.str = which(anno$strain == strains[ii])
    nreps[ii] = length(ind.str)

    #if there is more than one time when anno$strain is equal to strains[ii]
    #then take the average for the specific strain
    if(length(ind.str)>1){
      mat[,ii] = apply(expr[,ind.str],1,mean)
    } else {mat[,ii] = expr[,ind.str]}
  } # ii loop
  #Check if names are the same for each data frame
  names(mat) = strains
  rownames(mat) = rownames(expr)
  ann.out = cbind(strains,nreps) %>% as.data.frame(stringsAsFactors=F)
  return(list(expr=mat, anno=ann.out))
}

############################################################################

################################################################################
#########unique.genes.func
############################################

#' The function identofies probes that code for identical genes and takes the mean gene expression value across probes
#'
#' The first probe in the list is kept an contains reformatted, mean gene expression values
#' @param probe.vs.gene two column data frame containing the gene name assigned to each probe number
#' @param expr.mean data frame containing gene expression data for each strain of interests for a specific sex
#' @return a newly formatted data frame containing gene xpression data for each strain of interest
#' @return ind.test.4 with only probes for a unique gene (no duplicate genes)
#' @example
#' expr.mean.M.unique = unique.genes.func(probe.vs.gene = expr.formatted.gene.annotation, expr.mean = expr.mean.M) #Male strains
#'
unique.genes.func = function(probe.vs.gene = expr.formatted.gene.annotation, expr.mean = NULL) {

  #Mirror gene probes in between the key(=probe.vs.gene) and the expression data
  ind.test= probe.vs.gene[(probe.vs.gene$probe %in% rownames(expr.mean)),]

  #bind both data frames to have the probe and gene names
  ind.test1 = cbind(ind.test, expr.mean)

  #Groups rows that have the same gene names and takes the mean for every column of the duplicated rows
  library(dplyr)
  ind.test.2 = ind.test1 %>% group_by(gene) %>% mutate_each(list(mean), -(1)) %>% distinct

  #Eliminates all duplicated rows and keeps teh row and probe number of the first in the data frame
  ind.test.3 = ind.test.2[!duplicated(ind.test.2$gene),]

  #Reformatting
  #Deletes probe and gene column from the data frame
  #Uses probe as row names
  ind.test.4 = ind.test.3[,-c(1,2)]
  rownames(ind.test.4) = ind.test.3$probe

  return(ind.test.4)
}


##############################################################################
#########differences.func
#############################################################################

#' function to compute the Fold-Change in gene expression and phenotype between M and F
#'
#' (a substraction because the values are log transformed hence log(A)-log(B) = log (A/B))
#' Note that the function substracts the Male data from the Female. Hence a positive value means that the Female expression
#' or phenotype value is greater
#' !dPhno_log is unnecessary
#' @param expr.mean.F.unique mean Female expression data for every strain
#' @param expr.mean.M.unique mean Male expression data for every strain
#' @param pheno.mean.F mean Female phenotype data for every strain
#' @param pheno.mean.M mean Male phenotype data for every strain
#' @return a list containing 3 data frames
#' dExpr contains the FC of every gene expression for every strain
#' dPhno contains the FC of every phenotype for every strain
#' dPhno_log is the log of dPhno values
#' @example
#' differences.results=differences.func(expr.mean.F.unique = expr.mean.F.unique, expr.mean.M.unique = expr.mean.M.unique, pheno.mean.F=pheno.interesting.F,pheno.mean.M=pheno.interesting.M)

differences.func=function(expr.mean.F.unique=NULL,expr.mean.M.unique=NULL, pheno.mean.F=NULL, pheno.mean.M=NULL){

  # take differences, include logspace for phenotypes, F - M
  # column names are strains, rows are probes
  dExpr = expr.mean.F.unique - expr.mean.M.unique
  dPhno = pheno.mean.F - pheno.mean.M

  #unnecessary
  dPhno_log = log(pheno.mean.F) - log(pheno.mean.M)

  #Organize output as a list of 3 data frames
  out=list(dExpr=dExpr, dPhno=dPhno, dPhno_log=dPhno_log)
  return(out)
}

##################################################################################
##################cor.fn1
##################################################################################

#' Function for getting correlations for phenotype of interest FC and gene expression FC
#'
#' Need to compute dPhno for a specific phenotype prior
#' @param de = computed gene expression fold changes for every strain
#' @param dp = computed phenotype fold changes for every strain
#' @param thresh = correlation threshold
#' @param method is the correlation method. Here we use pearson correlation
#' @param p is the p value threshold, only data with p value smaller than p are computed
#' @return a x by 4 data frame with bicor values and p values for all x number of genes.
#' Specifies probe, phenotype, bicor, and p value
#' @example
#' cors1.NMR = cor.fn1(de=dExpr, dp=dPhno.NMR, thresh=0.4)
#'
cor.fn1 = function(de=NULL,dp=NULL,thresh=NULL,method="pearson",na.prop=0.4, p=1){

  start_time <- Sys.time()
  # empty vector
  res = c()

  #x is the the current row of the dExpr and y is the the current row of DPhno corresponds to a probe
  #go through each row
  for(ii in 1:nrow(de)){
    for(jj in 1:nrow(dp)){
      x = de[ii,] %>% data.matrix %>% as.numeric
      y = dp[jj,] %>% data.matrix %>% as.numeric
      # if there is a certain proportion of missing values, skip to next row
      if( length(which(is.na(x)==TRUE)) > na.prop*ncol(de) | length(which(is.na(y)==TRUE)) > na.prop*ncol(de) ){next}
      # compute correlation with pearson method, ignore na's on a casewise basis
      #rr = cor(x, y, method=method, use="complete.obs")
      #if correlation is na move on to next row
      rr=bicorAndPvalue(x,y,use="pairwise.complete.obs")
      if( is.na(rr$bicor) == TRUE ){next}
      #if correlation is higher than the threshold
      if(abs(rr$bicor) > thresh &(rr$p)<p){
        #add to the results vector a row with the row names corresponding to the Dexpr and Dphno being looked at
        # what is added would be probe, trait, correlation for the rows with significant correlation
        new = c(rownames(de)[ii], rownames(dp)[jj], rr$bicor, rr$p)
        res = rbind(res, new)

      }
    } # jj
  } # ii
  end_time <- Sys.time()

  #after it has gone through all of the rows, converts the results to a data frame
  res = as.data.frame(res, stringsAsFactors=F)
  names(res) = c("probe","pheno","bicor","p-value")
  rownames(res) = 1:nrow(res)
  res$bicor = data.matrix(res$bicor) %>% as.numeric

  return(res)

}

###############################################################
################################################################################################
#############find.gene.name.func
####################################################################################

#' Returns name of gene according to a probe number
#'
#' @param probe.number for which you want to know the gene name of
#' @return gene , the gene name
#' @example
#' find.gene.name.func("1451588_at")
#'
find.gene.name.func = function(probe = ""){

  gene = expr.formatted.gene.annotation$gene[which(expr.formatted.gene.annotation$probe == probe)]
  return(gene)
}


#########################################################################################################
####maxcor.gene.name.func
########################################################################################

#' Search for max absolute correlation and associated gene name
#'
#' function to search for probe and gene name of gene with max absolute correlation with a phenotype
#' @param cors1 data frame with correlation values for all strains. cors1."pheno" is specific for each phenotype
#' @return a row and a gene name
#' the row corresponds to the information of the gene with the max corr value
#' gene name is the name of the probe number found prior
#' @example
#' maxcor.gene.name.func(cors1 = cors1.NMR)
#'
maxcor.gene.name.func = function(cors1=NULL){

  #find max correlation value in data set
  ind = which(abs(cors1$bicor) == max(abs(cors1$bicor)))
  print(cors1[ind,])
  #find gene name associated to probe number
  gene = expr.formatted.gene.annotation$gene[which(expr.formatted.gene.annotation$probe == cors1[ind,]$probe)]
  return(gene)
}


################################################################################################
############select.FCcorG.func
###########################################################

#' Function to select genes with highest correlation
#'
#' Uses p.adjust to compute a false discovery rate for genes satisfying a correlation threshhold condition
#' and locate promising FCcorGenes
#' @param cors1 data frame containing the bicor and p-values for all genes for a specific phenotype
#' @param bicor.thresh minimum absolute bicor value |r| for consideration
#' @return cors1.FDR a data frame with the bicor and pvalues satisfying the bicor.thresh condition binded
#' with the FDR values for each gene that satisfied the condition
#' @example
#' results.NMR.select.FCcorG.func= select.FCcorG.func(cors1 = cors1.NMR, bicor.thresh = 0.4)
#'
select.FCcorG.func = function(cors1 = NULL, bicor.thresh = NULL){

  #select only genes from cors1 that satisfy the correlation requirement
  selected.bicor = subset(cors1 , abs(bicor)>bicor.thresh)

  #Isolated only p values
  pvals = selected.bicor$`p-value`

  #False discovery rate analysis (not needed)
  FDR = p.adjust(pvals, method = "BH", n = length(pvals))
  cors1.FDR = cbind(selected.bicor, FDR) %>% as.data.frame(stringsAsFactors=F)
  return(cors1.FDR)
}


#####################################################################
######### deg.analysis
#####################################################################

#' function for DEG analysis
#'
#' The function uses limma to calculate he values of differentially expressed genes in two data sets
#' @example
#' output = deg.analysis(dat_Null=expr.F.Null, dat_Alt=expr.M.Alt, null=null, alt=alt)
#'
deg.analysis = function(dat_Null=NULL,dat_Alt=NULL,null=NULL,alt=NULL){

  # implement linear model analysis with eBayes and BH adjustments
  subQall = rbind(dat_Null, dat_Alt)
  subQsex = c(rep(null,nrow(dat_Null)), rep(alt,nrow(dat_Alt)))
  design <- model.matrix(~0+subQsex)
  colnames(design) <- c(null,alt)
  contrast <- makeContrasts(Female - Male, levels = design)
  fit <- lmFit(t(subQall), design)
  fit <- contrasts.fit(fit, contrast) %>% eBayes
  output0 <- topTable(fit,number=ncol(subQall),adjust.method="BH")

  # convert log2FC to foldchange increase/decrease
  output = output0 %>% mutate(FC = 2^logFC)
  fc = rep(1,nrow(output))
  for (ii in 1:nrow(output)){
    if(output$FC[ii] > 1){fc[ii] = output$FC[ii]}
    if(output$FC[ii] < 1){fc[ii] = -1/output$FC[ii]}
  }
  output = output %>% mutate(ratioFC = fc) %>% mutate(absFC = abs(ratioFC))
  rownames(output) = rownames(output0)
  return(output)
}

####################################################################
###########select.deg.func
####################################################################

#' Function to select strong differentially expressed genes (DEG's)
#'
#' The function uses the output data from the prior executed deg.analysis function to select the genes with
#' significant DEG, above statistical significance. Threshholds such as the minimum Fold-CHange or a maximum
#' p-value can be set tot filter the data
#' @param data.deg the output from deg.analysis function
#' @param minlogFC.thresh let's use 1
#' @param minpvalue generally 0.05
#' @return selected.deg a data frame with the DEG genes that follow the criterias set
#' @example
#' selected.deg = select.deg.func(data.deg = output, minlogFC.thresh = 1, maxpvalue = 0.05)
#'
select.deg.func = function(data.deg = NULL, minlogFC.thresh = NULL, maxpvalue = NULL){

  probe.number = rownames(data.deg)
  selected.deg = cbind(probe.number, data.deg)
  selected.deg = selected.deg %>% filter(abs(logFC) > minlogFC.thresh, P.Value < maxpvalue) #Filter data with thresholds
  selected.deg = selected.deg[with(selected.deg, order(-logFC)), ]
  return(selected.deg)
}

#####################################################
################probe.to.gene.name.func
#####################################################

#' The function identifies probes that code for identical genes and takes the mean gene expression value across probes
#'
#' The first probe in the list is kept an contains reformatted, mean gene expression values
#' @param probe.vs.gene two column data frame containing the gene name assigned to each probe number
#' @param expr.mean data frame containing gene expression data for each strain of interests for a specific sex
#' @return a newly formatted data frame containing gene epression data for each strain of interest
#' @return ind.test.4 with only GENE NAMES for a unique gene (no duplicate genes). The gene names are trownames
#' @example
#' test.expr.mean.M.unique = probe.to.gene.name.func(probe.vs.gene = expr.formatted.gene.annotation, expr.mean = expr.mean.M) #Male strains
#'
probe.to.gene.name.func = function(probe.vs.gene = expr.formatted.gene.annotation, expr.mean = NULL) {

  #Mirror gene probes in between the key(=probe.vs.gene) and the expression data
  ind.test= probe.vs.gene[(probe.vs.gene$probe %in% rownames(expr.mean)),]

  #bind both data frames to have the probe and gene names
  ind.test1 = cbind(ind.test, expr.mean)

  #Groups rows that have the same gene names and takes the mean for every column of the duplicated rows
  library(dplyr)
  ind.test.2 = ind.test1 %>% group_by(gene) %>% mutate_each(list(mean), -(1)) %>% distinct

  #Eliminates all duplicated rows and keeps the row and probe number of the first appearance in the data frame
  ind.test.3 = ind.test.2[!duplicated(ind.test.2$gene),]

  #Reformatting
  #Deletes probe and gene column from the data frame
  #Uses probe as row names
  ind.test.4 = ind.test.3[,-c(1,2)]
  rownames(ind.test.4) = ind.test.3$gene

  return(ind.test.4)
}

##########################
