heatmap_RWHN <- function(rwhn_output, ylab, colours = c(low = "#e6e4f8", high = "#46009e")){
  
rwhn_flt <-  lapply(1:length(rwhn_output), function(i){
  rwhn_output[[i]] %>% 
    mutate(seed = i,
           rank = 1:nrow(.))
}) %>% do.call(rbind, .) %>% 
  group_by(name) %>% 
  mutate(rank_dif = abs(rank - mean(rank)),
        color = ifelse(rank_dif == 0, T, NA),
        V1 = signif(V1, digits = 2)) %>% 
  arrange(desc(V1)) %>% 
  #filter(rank_dif > (0.05 / max(rank_dif))) %>% 
  filter(rank_dif > 1) %>%                    # Filter terms that appear in the same position in all conditions
  ungroup() %>% 
  group_split(seed) %>% 
  lapply(., function(i){
    if(nrow(i > 0)){
      pct <- i[1,]$V1                             # Filter top 5% of terms ( with a messy loop!!)
      df <- i[1,]
      if(nrow(i) > 1){
        for(x in 2:nrow(i)){
          if(pct < 0.05){
            df <- rbind(df, i[x,])
            pct <- pct + i[x,]$V1
          }else{
            break
          }
        }
      }
     df$rank <- 1:nrow(df) 
    }else{
      df <- data.frame()
    }
    
    return(df)
  }) %>% 
  do.call(rbind,.) %>% 
  filter(!grepl("0\\.", name))

sighm <- ggplot(rwhn_flt, aes(x = as.factor(seed), y = name)) +
  geom_tile(aes(fill = rank), colour = "white") +
  scale_fill_gradient(breaks = seq(5, max(rwhn_flt$rank),5),
                      low = colours[1], high = colours[2]) +
  guides(color = FALSE,
         fill = guide_colourbar(title="Rank", reverse = T)) +
  xlab("Seed nodes") +
  ylab(ylab) +
  theme(legend.key.size = unit(.25, "cm"),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5), 
        axis.text.y = element_text(size = 5),
        axis.title = element_text(size = 6),
        panel.background = element_rect(fill = "black"), 
        panel.grid = element_blank(),
        legend.margin = margin(0,0,0,0, "cm")
  ) +
  scale_x_discrete(position = "top")  

return(sighm)
}
