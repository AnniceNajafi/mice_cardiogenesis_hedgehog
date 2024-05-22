#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Summer 2024
#' The following program is written to plot confidence intervals for different 
#' sliding window lengths



h <- 1
marker <- "Shh"

for(window_size in seq(10, 50, 5)){
  
  sliding_ligand.df <- marker_expression_SHF_OFT_traj(sliding_window_size = window_size, center_x, center_y, inner_radius, outer_radius, marker)
  
  
  bind_rows(sliding_ligand.df, .id = "source")->combined_df
  
  
  # Calculate the mean and confidence intervals
  summary_df <- combined_df %>%
    group_by(source) %>%
    summarise(
      mean_ligand = mean(ligand, na.rm = TRUE),
      lower_ci = mean_ligand - qt(0.975, df = n() - 1) * sd(ligand, na.rm = TRUE) / sqrt(n()),
      upper_ci = mean_ligand + qt(0.975, df = n() - 1) * sd(ligand, na.rm = TRUE) / sqrt(n())
    )
  
  # Plot using ggplot2
  hold.plt <- ggplot(summary_df, aes(x = as.numeric(source), y = mean_ligand, group = 1)) +
    geom_line(color = "#E1AFD1") +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), color="#AD88C6", alpha = 0.2) +
    labs(title = seq(10, 50, 5)[h],
         x = "Source",
         y = "Mean Ligand Value") +
    theme_classic()+
    #scale_color_manual(values = c("#E1AFD1", "#AD88C6"))+
    #geom_hline(yintercept =0)+
    theme(
      # Remove panel border
      panel.border = element_blank(),
      # Remove panel grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Remove panel background
      panel.background = element_blank(),
      #legend.position = "none",
      plot.title = element_text(hjust = 0.5, size=20),
      axis.text = element_blank(),
      text = element_text(size=18))+labs(x="Location from SHF-OFT", y=marker, tag=LETTERS[1:length(seq(10, 50, 5))][h])
  
  
      
      paste0("plt_", h)->name
      assign(name, hold.plt)
      h+1->h
}

(plt_1 + plt_2 + plt_3) / (plt_4 + plt_5 + plt_6) / (plt_7 + plt_8 + plt_9)















