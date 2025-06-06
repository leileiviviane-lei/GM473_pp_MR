

library(grid)
library(forestploter)

mydata=read.table("all_result.txt",header = T,sep = "\t")
mydata$` ` <- paste(rep(" ", 20), collapse = " ")
mydata$`OR (95% CI)` <- ifelse(is.na(mydata$or), "",sprintf("%.4f (%.4f - %.4f)",
                                                            mydata$or, mydata$or_lci95,
                                                            mydata$or_uci95))

forest(mydata[, c(2,4:5, 10:11)],
       est = mydata$or,
       lower = mydata$or_lci95,
       upper = mydata$or_uci95,
       sizes = 0.3,
       ci_column = 4,
       ref_line = 1,
       xlim = c(0.05, 3))
dev.copy2pdf(file = "forest_plot.pdf", width = 12, height = 6)
