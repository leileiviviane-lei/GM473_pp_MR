
library(circlize)
library(ComplexHeatmap)
library(tidyr)

df = read.csv("TableS5.csv", check.names = FALSE)
colnames(df)

# Data preprocessing
ivw = df[df$Method %in% c("Inverse variance weighted"), ]
ivw = ivw[, c("Exposure", "pval")]
df = df[, c("Exposure", "Method", "OR")]
df = pivot_wider(df, names_from = "Method", values_from = "OR")
df = merge(df, ivw, by = "Exposure")
rownames(df) = df$Exposure
df = df[, -1]

# Color mapping
color = list(
  col_1 = colorRamp2(c(0.47, 1.0, 1.91), c("#FF4500", "white", "#8A2BE2")),  # MR Egger
  col_2 = colorRamp2(c(0.69, 1.0, 1.66), c("#00CED1", "white", "#FF1493")),  # Weighted median
  col_3 = colorRamp2(c(0.73, 1.0, 1.65), c("#7FFF00", "white", "#DC143C")),  # Inverse variance weighted
  col_4 = colorRamp2(c(0.51, 1.0, 2.53), c("#FF8C00", "white", "#1E90FF")),  # Simple mode
  col_5 = colorRamp2(c(0.62, 1.0, 1.83), c("#FFD700", "white", "#8B008B")),  # Weighted mode
  col_6 = colorRamp2(c(0.001, 0.025, 0.05), c("#9400D3", "white", "#00FFFF")) # p-value
)

# Prepare inner circle data (columns 1:5)
aa = as.matrix(df[, 1:5, drop = FALSE])

circos.clear()
circos.par(gap.degree = 1)
circos.par$start.degree = 20
circos.par(track.margin = c(0.001, 0.001))

# Plot inner circle (5 tracks with different colors)
for (i in 1:ncol(aa)) {
  circos.heatmap(as.matrix(aa[, i, drop = FALSE]),
                 col = color[[i]],
                 cluster = FALSE,
                 track.height = 0.2 / ncol(aa),
                 cell.lwd = 0.5,
                 cell.border = "white")
}

# Prepare outer circle data (6th column, p-value)
bb <- as.matrix(df[, 6, drop = FALSE])
mode(bb) <- "numeric"

# Plot outer circle (p-value track)
circos.heatmap(bb,
               col = color[[6]],
               cluster = FALSE,
               track.height = 0.1,
               cell.lwd = 0.5,
               cell.border = "white",
               rownames.side = "outside",
               rownames.cex = 0.25)

# Legend setup
library(ComplexHeatmap)
library(grid)

lgd1 <- Legend(title = "MR Egger", 
               border = "black", 
               grid_height = unit(3, "mm"), 
               legend_width = unit(30, "mm"), 
               at = c(0.47, 1.0, 1.91),  
               title_position = "topcenter", 
               col_fun = color[[1]], 
               direction = "horizontal",
               title_gp = gpar(fontsize = 10, fontface = "bold"))

lgd2 <- Legend(title = "Weighted median", 
               border = "black", 
               grid_height = unit(3, "mm"), 
               legend_width = unit(30, "mm"), 
               at = c(0.69, 1.0, 1.66),  
               title_position = "topcenter", 
               col_fun = color[[2]], 
               direction = "horizontal",
               title_gp = gpar(fontsize = 10, fontface = "bold"))

lgd3 <- Legend(title = "Inverse variance weighted", 
               border = "black", 
               grid_height = unit(3, "mm"), 
               legend_width = unit(30, "mm"), 
               at = c(0.73, 1.0, 1.65),  
               title_position = "topcenter", 
               col_fun = color[[3]], 
               direction = "horizontal",
               title_gp = gpar(fontsize = 10, fontface = "bold"))

lgd4 <- Legend(title = "Simple mode", 
               border = "black", 
               grid_height = unit(3, "mm"), 
               legend_width = unit(30, "mm"), 
               at = c(0.51, 1.0, 2.53),  
               title_position = "topcenter", 
               col_fun = color[[4]], 
               direction = "horizontal",
               title_gp = gpar(fontsize = 10, fontface = "bold"))

lgd5 <- Legend(title = "Weighted mode", 
               border = "black", 
               grid_height = unit(3, "mm"), 
               legend_width = unit(30, "mm"), 
               at = c(0.62, 1.0, 1.83),  
               title_position = "topcenter", 
               col_fun = color[[5]], 
               direction = "horizontal",
               title_gp = gpar(fontsize = 10, fontface = "bold"))

lgd6 <- Legend(title = "p-value", 
               border = "black", 
               grid_height = unit(3, "mm"), 
               legend_width = unit(30, "mm"), 
               at = c(0.001, 0.025, 0.05),  
               title_position = "topcenter", 
               col_fun = color[[6]], 
               direction = "horizontal",
               title_gp = gpar(fontsize = 10, fontface = "bold"))

# Combine all legends
pd <- packLegend(lgd1, lgd2, lgd3, lgd4, lgd5, lgd6, row_gap = unit(2, "mm"))

# Export legends to PDF
pdf("Legend_Plot.pdf", width = 6, height = 4)
draw(pd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
dev.off()

# Export legends to TIFF
tiff("Legend_Plot.tiff", width = 1800, height = 1200, res = 300)
draw(pd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
dev.off()
