# Function to plot binary adjacency matrix as a heatmap.

# @param G.mat       binary adjacency matrix (square matrix with 0/1 entries)
# @param title       character; title for the plot
# @param output_path character; path to save the plot
# @param plot_name   character; name of the plot file
# @return           saves a heatmap plot of the adjacency matrix to the specified path
plot_binary_adj_matrix <- function(G.mat, title, output_path, plot_name){
  adj_df <-  melt(G.mat)
  colnames(adj_df) <- c("Row", "Col", "Value")
  adj_df$Value <- as.numeric(adj_df$Value)
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "grey92", linewidth = 0.15) +
    coord_fixed() +
    scale_y_reverse() +
    labs(
      title = str_wrap(title, width = 42),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(
        size = 20,
        hjust = 0.5,
        vjust = -1.5,
        lineheight = 0.95 
      ),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.35, "cm"),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    aes(fill = factor(Value, levels = c(0, 1))) + scale_fill_manual(
      values = c( "0" = "white","1" = "black"),
      labels = c("0" = "Absent", "1" = "Present"),
      name = "Edge"
    )  
  ggsave(plot, filename=paste(output_path, plot_name, sep=""), width=8, height=8, dpi=300)
}
