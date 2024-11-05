pairwise_correlations <- function(df, group_var) {
  
  
  df %>%
    # Group by the group_var (can be any variable like 'group', 'ID', etc.)
    group_by({{ group_var }}) %>%
    # Nest data by group to iterate over each group separately
    nest() %>%
    # Apply pairwise correlation within each group
    mutate(cor_results = map(data, ~{
      data <- .x
      # Get variable names
      var_names <- colnames(data)
      # Get all pairwise combinations of variables
      combs <- combn(var_names, 2, simplify = FALSE)
      
      # Perform pairwise correlations for all combinations
      map_dfr(combs, function(pair) {
        cor_test <- cor.test(data[[pair[1]]], data[[pair[2]]], use = "complete.obs")
        tibble(
          var1 = pair[1],
          var2 = pair[2],
          correlation = cor_test$estimate,
          p_value = cor_test$p.value
        )
      })
    })) %>%
    # Unnest the correlation results to get a flat table
    unnest(cor_results)
}

pw.result <- pairwise_correlations(tools_wide, group_var = Samples)

library(RColorBrewer)
library(hrbrthemes)
library(patchwork)

hm_list = list()
rng = range(pw.result$correlation)

# no x y labs
# no grid lines
# vertical x axis lab
# smaller title
# print rounded number in box

for(s in unique(pw.result$Samples)){
  
  hm_list[[s]] = 
    pw.result %>%
    filter(Samples == s) %>%
    mutate(CorrPrint = round(correlation, 2)) %>%
    ggplot(aes(var1, var2, fill= correlation)) + 
    geom_tile(color = "white") +
    geom_text(aes(label = CorrPrint), 
              color = "white",
              size = 3) +
    coord_equal() +
    scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd"),
                         limits=c(floor(rng[1]), ceiling(rng[2]))) +
    ggtitle(s) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "bold")) 
}

wrap_plots(hm_list) + guide_area() + plot_layout(guides = 'collect')
ggsave(paste0(results.dir, "/correlation_heatmap.png"),
       width = 15,
       height = 12)