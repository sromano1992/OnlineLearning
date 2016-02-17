plotTable<-function (table,fileName){  
  tg <- tableGrob(table, rows = seq_len(nrow(table))) 
  
  fullheight <- convertHeight(sum(tg$heights), "cm", valueOnly = TRUE)
  margin <- convertHeight(unit(0.51,"in"), "cm", valueOnly = TRUE)
  a4height <- 29.7 - margin
  nrows <- nrow(tg)
  npages <- ceiling(fullheight / a4height)
  
  heights <- convertHeight(tg$heights, "cm", valueOnly = TRUE) 
  rows <- cut(cumsum(heights), include.lowest = FALSE,
              breaks = c(0, cumsum(rep(a4height, npages))))
  
  groups <- split(seq_len(nrows), rows)
  
  gl <- lapply(groups, function(id) tg[id,])
  
  pdf(fileName, paper = "special", width = 14.8, height = 21)  #A5
  for(page in seq_len(npages)){
    grid.newpage()
#     grid.rect(width=unit(21,"cm") - unit(0.51,"in"),
#               height=unit(29.7,"cm")- unit(0.51,"in"))
    grid.draw(gl[[page]])
  }
  ## alternative to explicit loop:
  ## print(marrangeGrob(grobs=gl, ncol=1, nrow=1, top=NULL))
  dev.off()
}