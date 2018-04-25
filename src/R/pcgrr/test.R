library(DT)
options(DT.options = list(pageLength = 5))
df = as.data.frame(cbind(matrix(round(rnorm(50), 3), 10), sample(0:1, 10, TRUE)))
# style V6 based on values of V6
datatable(df, extensions='Responsive') %>% formatStyle(
  'V6',
  backgroundColor = styleEqual(c(0, 1), c('gray', 'yellow'))
)
