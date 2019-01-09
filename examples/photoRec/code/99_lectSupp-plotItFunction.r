## makes stripplot for a handful of probesets ... used in all my lecture support
## scripts that show what different kinds of DEA "hits" look like

# stripplotIt <- function(myData, ...) {
#     stripplot(gExp ~ devStage | gene, myData,
#               group = gType, jitter.data = TRUE,
#               auto.key = TRUE, type = c('p', 'a'), grid = TRUE, ...)
# }

stripplotIt <- function(myData, yFree = FALSE, ...) {
  require(lattice)
  thePlot <- stripplot(gExp ~ devStage | gene, myData,
                       type = c('p', 'a'), grid = TRUE,
                       group = gType, auto.key = TRUE,
                       jitter.data = TRUE, ...)
  if(yFree) {
    thePlot <- update(thePlot,
                      scales = list(y = list(relation = "free")))
  }
  print(thePlot)
}

## similar but using xyplot()
xyplotIt <- function(myData, yFree = FALSE, myType = c('p', 'a'), ...) {
  require(lattice)
  thePlot <- xyplot(gExp ~ age | gene, myData,
                    type = myType, grid = TRUE,
                    group = gType, auto.key = TRUE,
                    jitter.data = TRUE, ...)
  if(yFree) {
    thePlot <- update(thePlot,
                      scales = list(y = list(relation = "free")))
  }
  print(thePlot)
}