plot_grid(pf1 + ggtitle("Single-cell run 1"),pf2 + ggtitle("Single-cell run 2"),labels=c("A","B"))
ggsave("fig_final//s5_marker_genes.png",bg="white", width=16,height=12,dpi=300)
ggsave("fig_final/svg/s5_marker_genes.svg",bg="white", width=16,height=12,dpi=300)
