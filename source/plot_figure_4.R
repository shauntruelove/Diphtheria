
# Source Functions to make Figure 4
source('source/figure_4_functions.R')



# PLOT ALL TOGETHER -------------------------------------------------------

A <- plot_combined_VcR0(delay=1, alpha=0, alpha_density_R0=0.5)
B <- plot_combined_VcR0(delay=1, alpha=0.25, alpha_density_R0=0.5)
C <- plot_combined_VcR0(delay=2, alpha=0.5, alpha_density_R0=0.5)
D <- plot_combined_VcR0(delay=5, alpha=.9, alpha_density_R0=0.5)

plot_grid(A, B, C, D, nrow=2, ncol=2, labels=c('A','B','C','D'))


# Save PDF of figure 4

png('results/figure4.png', width = 1225, height=945)
plot_grid(A, B, C, D, nrow=2, ncol=2, labels=c('A','B','C','D'))
dev.off()

