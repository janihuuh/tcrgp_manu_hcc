


list.of.packages <- c("ggplot2", "dplyr", "plyr", "here", "data.table", "RColorBrewer", "Seurat", "ComplexHeatmap", "circlize")
new.packages     <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  
  message(paste("Installing missing packages", new.packages))
  install.packages(new.packages)
  
}

message("All packages installed")
