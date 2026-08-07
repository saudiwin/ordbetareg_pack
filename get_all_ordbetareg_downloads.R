library(dplyr)

get_d <- lapply(seq.Date(from=as.Date("2022-07-01"),to=as.Date("2025-09-23"),by=1), function(i) {
  
  
  c1 <- try(packageRank::packageRank("ordbetareg",date=i))
  
  if('try-error' %in% class(c1)) return(tibble())
  
  c1$package.data
  
  
}) %>% bind_rows
