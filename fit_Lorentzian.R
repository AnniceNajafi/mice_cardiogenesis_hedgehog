#' Authors@R: person("Annice", "Najafi", email = "annicenajafi@tamu.edu")
#' Department of Biomedical Engineering, Texas A&M University, College Station, TX
#' Summer 2024
#' The following program is written to Fit Lorentzian to the ligand data

#Load library 
#Note: nls doesn't work and gives errors, try minpack.lm
library(minpack.lm)

mean_ligand_df <- combined_df %>% group_by(source) %>% summarise(ligand = mean(ligand, na.rm = TRUE))

#Convert to numeric
as.numeric(mean_ligand_df$ligand)->mean_ligand_df$ligand
as.numeric(mean_ligand_df$source)->mean_ligand_df$source


# Determine initial estimates
initial_a <- max(mean_ligand_df$ligand)
initial_c <- mean_ligand_df$source[which.max(mean_ligand_df$ligand)]
initial_b <- diff(range(mean_ligand_df$source)) / 2  # A rough estimate

my_fit1 <- nlsLM(
  ligand ~ a * b^2 / ((source - c)^2 + b^2),
  data = mean_ligand_df,
  start = list(a = initial_a, c = initial_c, b = initial_b)
)


summary(my_fit1)


# Extract the estimated parameters
params <- coef(my_fit1)
a_est <- params["a"]
b_est <- params["b"]
c_est <- params["c"]


fitted_ligand <- a_est * b_est^2 / ((mean_ligand_df$source - c_est)^2 + b_est^2)
fitted_curve_df <- data.frame(source = mean_ligand_df$source , ligand = fitted_ligand)


ggplot() +
  geom_point(data = mean_ligand_df, aes(x = source, y = ligand), color = "#4D869C", shape=8, stroke=4) +
  geom_line(data = fitted_curve_df, aes(x = source, y = ligand), color = "#32012F", size=2) +
  labs(title = "Lorentzian Fit to Data",
       x = "Source",
       y = "Shh") +
  theme_classic()+
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
    text = element_text(size=18))+labs(x="Location from SHF-OFT", y="Shh", tag="B")




####Try it on the entire dataset


#Convert to numeric
as.numeric(combined_df$ligand)->combined_df$ligand
as.numeric(combined_df$source)->combined_df$source




# Determine initial estimates
initial_a <- max(combined_df$ligand)
initial_c <- combined_df$source[which.max(combined_df$ligand)]
initial_b <- diff(range(combined_df$source)) / 2  # A rough estimate

my_fit1 <- nlsLM(
  ligand ~ a * b^2 / ((source - c)^2 + b^2),
  data = combined_df,
  start = list(a = initial_a, c = initial_c, b = initial_b)
)


summary(my_fit1)


# Extract the estimated parameters
params <- coef(my_fit1)
a_est <- params["a"]
b_est <- params["b"]
c_est <- params["c"]


fitted_ligand <- a_est * b_est^2 / ((combined_df$source - c_est)^2 + b_est^2)
fitted_curve_df <- data.frame(source = combined_df$source , ligand = fitted_ligand)


ggplot() +
  geom_point(data = combined_df, aes(x = source, y = ligand), color = "#4D869C", shape=8, stroke=4) +
  geom_line(data = fitted_curve_df, aes(x = source, y = ligand), color = "#32012F", size=2) +
  labs(title = "Lorentzian Fit to Data",
       x = "Source",
       y = "Shh") +
  theme_classic()+
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
    text = element_text(size=18))+labs(x="Location from SHF-OFT", y="Shh", tag="B")

