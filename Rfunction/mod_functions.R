### This file contains modified functions from the MAGeCKFlute package. 
### The modification and reason for doing so is added as comment above each function.


# add bar.position arugment to default barplot functions to allow reordering of the bars (baseline T0 samples first).
MapRatesView <- function(countSummary,
                         Label = "Label",
                         Reads = "Reads",
                         Mapped = "Mapped",
                         bar.position = NULL,
                         filename = NULL,
                         width = 5, height = 4,
                         ...){
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package \"scales\" is required. Please install it.", call. = FALSE)
  }
  gg = data.frame(Label=rep(countSummary[, Label], 2),
                  read=rep(countSummary[, Reads],2),
                  count=c(countSummary[, Mapped],
                          countSummary[, Reads]-countSummary[, Mapped]),
                  category=factor(rep(c("mapped", "unmapped"),
                                      c(nrow(countSummary), nrow(countSummary))),
                                  levels = c("unmapped", "mapped")))
  gg$percent = paste0(round(gg$count*100/gg$read, 1), "%")
  gg$pos = ceiling(gg$count/2)
  gg$pos[gg$category=="unmapped"] = ceiling(gg$pos[gg$category=="unmapped"] +
                                              gg$pos[gg$category=="mapped"]*2)
  
  if(is.null(bar.position)){
    bar.position = Label
  }  # added
  fill = c("#9BC7E9", "#1C6DAB")
  
  p <- ggplot(gg)
  p = p + geom_bar(aes_string(y = "count", x = "Label", fill = "category"),
                   stat="identity", width=0.8, alpha=0.9)
  p = p + geom_text(aes_string(x = "Label", y = "pos", label = "percent"), size=4)
  p = p + labs(x=NULL, y="Reads", title="Mapping ratio")
  p = p + scale_y_continuous(expand = c(0,0))
  p = p + scale_x_discrete(limits = bar.position)  # added
  p = p + scale_fill_manual(values=fill)
  p = p + theme(legend.title = element_blank())
  p = p + theme(text = element_text(colour="black",size = 14),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text.x = element_text(angle = 45, hjust=1, vjust = 1,
                                           colour="gray10", face="plain"),
                axis.text.y= element_text(colour="gray10", face="plain"))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}




IdentBarView <- function(gg, x = "x", y = "y", fill = c("#CF3C2B", "#394E80"),
                         bar.position = NULL,
                         main = NULL, xlab = NULL, ylab = NULL,
                         filename = NULL, width = 5, height = 4, ...){
  gg$x = gg[, x]
  gg$y = gg[, y]
  
  if(is.null(bar.position)){
    bar.position = x
  }  # added
  
  p <- ggplot(gg)
  p = p + geom_bar(aes(x, y), stat="identity", width=0.6, fill = fill[1], alpha=0.8)
  p = p + labs(x=xlab, y=ylab, title=main)
  p = p + scale_y_continuous(expand = c(0,0))
  p = p + scale_x_discrete(limits = bar.position)  # added
  p = p + theme(text = element_text(colour="black",size = 14),
                plot.title = element_text(hjust = 0.5, size=18),
                axis.text = element_text(colour="gray10"),
                axis.text.x=element_text(angle = 45, hjust=1, vjust = 1))
  p = p + theme(axis.line = element_line(size=0.5, colour = "black"),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank())
  
  if(!is.null(filename)){
    ggsave(plot=p, filename=filename, units = "in", width=width, height=height, ...)
  }
  return(p)
}




CutoffCalling=function(d, scale=1){
  param=1
  if(is.logical(scale) & scale){
    param = round(length(d) / 20000, digits = 1)
  }else if(is.numeric(scale)){param = scale}
  
  Control_mean=0
  sorted_beta=sort(abs(d))
  temp=quantile(sorted_beta,0.68)
  temp_2=qnorm(0.84)
  cutoff=round(temp/temp_2,digits = 3)
  names(cutoff)=NULL
  cutoff=cutoff*param
  return(cutoff)
}