

tmp <- read.csv2("C:/Users/jerem/Downloads/MSFT_2018_2022_30min.txt",sep=',')
tmp <- tmp[,c("X.DATE.","X.TIME.","X.CLOSE.")]
colnames(tmp) <- c('date','time','close')
tmp$time <- as.character(tmp$time)
tmp$time <- ifelse(nchar(tmp$time)<=5,
                   paste0("0",tmp$time),
                   tmp$time)
tmp[["dateTime"]] <- as.Date(paste(tmp$date,tmp$time),"%Y%m%d %H%M%S")


tmp_fd <-
  tmp %>%
  select('date','time','close') %>%
  pivot_wider(id_cols = 'time', id_expand = TRUE,
              names_from = 'date', values_from = 'close') %>%
  filter(!row_number() %in% nrow(.)) %>%
  select(-c(time)) %>%
  mutate(across(1:ncol(.), as.numeric)) %>%
  # mutate(across(`20180102`:`20221230`, as.numeric)) %>%
  as.data.frame() %>%
  linear_imputatation()

tmp_cidr <- compute_cidr(tmp_fd)

plot_fd(tmp_fd)
plot_fd(tmp_cidr)

saveRDS(tmp_fd,"C:/Users/jerem/Downloads/msft_cidr_30.rds")
saveRDS(tmp_cidr,"C:/Users/jerem/Downloads/msft_base_30.rds")
