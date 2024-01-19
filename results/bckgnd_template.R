library(stringr)
library(lubridate)
library(job)


tstamp <- now() %>% toString() %>% str_replace(' ','___')
tstamp