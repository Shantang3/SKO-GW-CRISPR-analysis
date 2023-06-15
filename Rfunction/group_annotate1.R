group_annotate <- function(filename, norm.method="loess", 
                          samples, use.diff=TRUE,sd){
  
  require(dplyr)
  
  set.seed(42)
  
  mle <- read.delim(filename, header = T, sep = "\t", stringsAsFactors = F)
  mle <- mle %>% dplyr::select(-contains(".beta"), -contains(".wald"), -contains(".z"))
  
  beta <- MAGeCKFlute::ReadBeta(filename)
  beta_norm <- MAGeCKFlute::NormalizeBeta(beta, samples = samples, method = norm.method)
  
  dd <- inner_join(beta_norm, mle)
  dd$diff <- dd[,samples[2]]-dd[,samples[1]]
  
  
  # annotate groups
  x_cutoff=CutoffCalling(dd[,samples[1]],sd)
  x_cutoff=sort(c(-x_cutoff, x_cutoff))
  y_cutoff=CutoffCalling(dd[,samples[2]],sd)
  y_cutoff=sort(c(-y_cutoff, y_cutoff))
  intercept=CutoffCalling(dd$diff,sd)
  
  if(use.diff){
    idx0 <- dd$diff<intercept & dd$diff>(-intercept)
  }else idx0=0
  
  dd$group="Others"
  idx1=dd[,samples[1]] < x_cutoff[1]
  idx2=dd[,samples[1]] > x_cutoff[2]
  idx3=dd[,samples[2]] < y_cutoff[1]
  idx4=dd[,samples[2]] > y_cutoff[2]
  
  dd$group[idx1&idx3]="bottomleft"
  dd$group[idx1&idx4]="topleft"
  dd$group[idx2&idx3]="bottomright"
  dd$group[idx2&idx4]="topright"
  dd$group[!idx1&!idx2&idx3]="bottomcenter"
  dd$group[!idx1&!idx2&idx4]="topcenter"
  dd$group[!idx3&!idx4&idx1]="midleft"
  dd$group[!idx3&!idx4&idx2]="midright"
  
  grp.interest <- c("midleft", "bottomcenter", "midright", "topcenter")
  dd$group[!dd$group%in%grp.interest]="Others"
  dd$group[idx0]="Others"
  
  return(dd)
}








