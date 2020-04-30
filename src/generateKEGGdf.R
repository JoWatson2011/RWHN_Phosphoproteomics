library(dplyr)

kegg <- read.csv("Downloads/KEGG_2019_Human.txt", col.names = c("tmp"), stringsAsFactors = F)

rmv <- kegg %>% filter(!grepl("\t", tmp) | grepl("^ ", tmp))
kegg <- kegg %>% 
  filter(!grepl("\t", tmp) | grepl("^ ", tmp)) %>% 
  mutate(gr = rep(1:(nrow(rmv)/2), 2)[order(rep(1:(nrow(rmv)/2), 2))]) %>% 
  group_by(gr) %>% 
  summarise(tmp = paste(tmp, collapse = "")) %>%
  select(tmp) %>% 
  rbind(kegg) %>% 
  filter(!(tmp %in% rmv$tmp)) %>% 
  separate(tmp, c("pw", "genes"), sep = "\t\t") %>% 
  separate_rows(genes, sep = "\t") %>% 
  group_by(pw) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  #THIS IS NOT A P VALUE... 
  # More like a observation value?!
  mutate(p = n/sum(n)) %>% 
  select(-n)

readr::write_tsv(kegg, "kegg.txt", col_names = F)
