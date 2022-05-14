evoFreq_fixColor_clone <- function(df) {
  
  # Set color scale
  palette <- c(palette.colors(palette = "R4"), "#FFC0CB", "#006400")
  
  for (i in 1:length(palette)) {
    df[df$clone_id==(i+min(df$clone_id)-1),"plot_color"] <- palette[i]
  }
  return(df)
}



evoFreq_fixColor_clus <- function(df) {
  
  # Set color scale
  palette <- c("#1F78B4", "#33A02C", "#A6CEE3", "#E31A1C", "#6A3D9A", "#B2DF8A", 
               "#FB9A99", "#000000", "#FF7F00",  "#B15928", "#FFFF99", "#CAB2D6", "#FDBF6F")
  
  for (i in 1:length(palette)) {
    df[df$clone_id==(i-1),"plot_color"] <- palette[i]
  }
  return(df)
}
