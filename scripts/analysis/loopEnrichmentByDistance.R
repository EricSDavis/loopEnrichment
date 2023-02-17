## Script to understand how loop
## enrichment changes with distance

## Load required packages
library(mariner)
library(data.table)
library(InteractionSet)

## Define file paths to Hi-C and loop files
hicFile <- "data/raw/HPOW_data/hic/rao_all/HPOW_tot_inter.hic"
loopFile <- "data/raw/HPOW_data/hic/rao_all/HPOW_tot_5kbLoops.txt"

## Define shared params
binSize <- 10e03
buffer <- 10

## Read in loops
loops <- 
  loopFile |>
  fread() |>
  binPairs(binSize)

## Check that all data is in the upper triangular
all(start1(loops) < start2(loops))

## Extract and calculate loop enrichment ---------------------------------------

## Extract matrices around each loop
mats <- 
  pixelsToMatrices(x=loops, buffer=buffer) |>
  pullHicMatrices(
    files=hicFile,
    binSize=binSize,
    half="upper"
  )

## Don't need to filter out loops b/c
## Selecting TopLeft and BottomRight
## as background is parallel to diagonal, so
## we can use short loops also.
# ## Filter out loops that cross the diagonal
# keep <- apply(counts(mats), 3, \(x) {
#   !anyNA(x[,,1])
# })
# mats <- mats[keep]
# 
# ## What perc of loops are excluded
# table(pairdist(loops) < 200e3)

## Define center pixels
cp <- counts(mats)[buffer+1, buffer+1,,]

## Define background
idx1 <- mariner:::.selectTopLeft(buffer=buffer, 4, invert=FALSE)
idx2 <- mariner:::.selectBottomRight(buffer=buffer, 4, invert=FALSE)
idx <- c(idx1, idx2)
bg <- apply(counts(mats), c(3,4), \(x) median(x[idx] + 1))


## Add loop enrichment score and pairdistance
interactions(mats)$score <- cp / bg[,1]
interactions(mats)$size <- pairdist(mats)

## Function to calculate loop enrichment (from mariner)
# enrich <- calcLoopEnrichment(
#   x=loops,
#   files=hicFile,
#   mhDist=c(4,5,6)
# )

## Visualize trend -------------------------------------------------------------

## Calculate rolling enrichment at different window sizes
## and correcting them.
ns <- round(seq(1, 200, length.out = 9))
par(mfrow=c(3,3))
for (n in ns) {
  .correctEnrich(x=interactions(mats),
                 scores=cp/bg[,1],
                 k=n,
                 nknots=10,
                 plot=TRUE)
}
par(mfrow=c(1,1))


