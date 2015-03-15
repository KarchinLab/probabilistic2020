##############################
# Calculates that statistical significance of frameshift mutations
# in genes by a likelihood ratio test compared to a non-coding
# background rate.
# 
# Parameters
# background : file with non-coding frameshift counts
# counts : file with frameshift counts in genes
# output : path to output results
###############################
suppressPackageStartupMessages(library(bbmle))
library(reshape2)
suppressWarnings(library(emdbook))

if ("getopt" %in% rownames(installed.packages())){
  # get command line arguments
  library(getopt)
  spec <- matrix(c(
    'background', 'b', 1, 'character',
    'counts', 'c', 1, 'character',
    'output', 'o', 1, 'character'
    ), byrow=TRUE, ncol=4)
  opt = getopt(spec)
  # print out help msg
  if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
  }
} else {
  opt <- list(ARGS=NULL)
}

frameshift.lrt <- function(bg, fs.counts){
  # real data
  df <- read.delim(bg, sep="\t", row.names=1)
  
  # get the number of nucleotides
  n <- df[1,"N"]
  
  # get the columns for frameshift counts
  indel.cols <- colnames(df)[grepl('^X', colnames(df))]
  tmp.df <- as.data.frame(df[,indel.cols])
  colnames(tmp.df) <- indel.cols

  # get the method of moment estimates for initialization of mle solving
  prob.moment <- colMeans(tmp.df) / n
  sample.variance <- apply(tmp.df, 2, var)
  theta.moment <- (sample.variance - n*prob.moment*(1-prob.moment)) / (n^2*prob.moment*(1-prob.moment)-sample.variance)
  
  # re-format data frame
  bg.indel.df <- melt(df, 
                      measure.vars=indel.cols,
                      variable.name='indel.len',
                      value.name='num.indels')
  
  # figure out mle of background
  # NOTE: prob.mle and theta.mle are considered global vars
  # sint this is needed for the mle2 formula
  theta.mle <<- theta.moment  # filler for mle vals until fitted
  prob.mle <<- prob.moment  # filler for mle vals until fitted
  for (i in 1:length(indel.cols)){
    # fit beta-binomial
    bgFit <- mle2(num.indels~dbetabinom(prob, size=N, theta=1/mytheta),   #theta),
                  parameters=list(prob~1,mytheta~1),
                  data=bg.indel.df[bg.indel.df["indel.len"]==indel.cols[i],],
                  start=list(prob=prob.moment[i], mytheta=theta.moment[i]),
                  method="L-BFGS-B",
                  control=list(parscale=unlist(list(prob=prob.moment[i], mytheta=theta.moment[i]))),
                  lower=c(prob=rep(3e-16, by=length(prob.moment[i])), mytheta=rep(3e-16, by=length(prob.moment[i]))))
  
    # record parameters
    theta.mle[i] <- coef(bgFit)[["mytheta"]]
    prob.mle[i] <- coef(bgFit)[["prob"]]
  }
  
  # estimate relative rate
  count.df <- read.delim(fs.counts, sep="\t", row.names=1)
  count.df <- count.df[count.df["total"]>0,]
  # check if there is data
  if (nrow(count.df)==0){
      return(NULL)
  }
  count.df["ratio.moment"] <- count.df["total"] / (count.df["bases.at.risk"]*sum(prob.mle))
  mycols <- c(indel.cols, "bases.at.risk", "ratio.moment")
  mycount.df <- count.df[,mycols]
  mycount.df["ratio.mle"] <- 1.0
  mycount.df["frameshift.p.value"] <- 1.0
  failedConvergence <- rep(0, by=nrow(mycount.df))
  for (i in 1:nrow(mycount.df)){
    # get indel counts for gene
    tmp.df <- mycount.df[i,]
    gene.indel.df <- melt(tmp.df, 
                          measure.vars=indel.cols,
                          variable.name='indel.len',
                          value.name='num.indels')
    
    # get MLE
    tryCatch(alt <- mle2(num.indels~dbetabinom(prob=ratio*prob.mle, size=bases.at.risk, theta=1/theta.mle),
                         parameters=list(ratio~1),
                         data=gene.indel.df,
                         start=list(ratio=gene.indel.df[1,"ratio.moment"]),
                         method="L-BFGS-B",
                         control=list(parscale=unlist(list(ratio=gene.indel.df[1,"ratio.moment"]))),
                         lower=c(ratio=3e-16)), 
             error=function(e) print(paste('ERROR:', row.names(tmp.df))))
    
    # get log likelihood ratio test statistic
    converge.status <- slot(alt, "details")$converge
    if (!converge.status){
      logLikAlt <- as.numeric(logLik(alt))
    }else{
      failedConvergence[i] <- converge.status
      print(paste(row.names(tmp.df), 'failed to converge. Reverting to method of moments estimate.'))
      logLikAlt <- sum(dbetabinom(gene.indel.df[,"num.indels"], gene.indel.df[1,"ratio.moment"]*prob.mle, gene.indel.df[,"bases.at.risk"], 1/theta.mle, log=T))
    } 
    logLikNull <- sum(dbetabinom(gene.indel.df[,"num.indels"], prob.mle, gene.indel.df[,"bases.at.risk"], 1/theta.mle, log=T))
    lrtStatistics <- 2*(logLikAlt-logLikNull)
    
    # get p-value
    pval <- 1-pchisq(lrtStatistics, 1)
    mycount.df[i,"ratio.mle"] <- coef(alt)[[1]]
    mycount.df[i,"frameshift.p.value"] <- pval
  }
  output <- mycount.df[order(mycount.df['frameshift.p.value']),]
  output['frameshift.q.value'] <- p.adjust(output[,'frameshift.p.value'])
  return(output)
}

if (length(opt)==4){
  # compute results
  bg <- opt$background
  fs.counts <- opt$counts
  print('Computing Likelihood Ratio Test . . .')
  output <- frameshift.lrt(bg, fs.counts)
  print('Finished Likelihood Ratio Test')
  output.path <- opt$output
  
  # re-order column output
  if (!is.null(output)){
    output$gene.name <- rownames(output)
    rownames(output) <- NULL
    colName <- colnames(output)
    output <- output[c(colName[ncol(output)], colName[1:(ncol(output)-1)])]
    
    # write output
    write.table(output, output.path, sep='\t', quote=F, row.names=F)
  } else {
    print('There were no frameshifts observed.')
  }
}
