path <- 'C:/Users/jerem/Downloads/'

x <- eseis::read_mseed(
    paste0(path, 'fdsnws-dataselect_2024-08-01t05_44_23z.mseed')
  )
plot(x)
