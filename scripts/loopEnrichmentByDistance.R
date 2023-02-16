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

## Visualize trend -------------------------------------------------------------

## Data table of loop size and score
dat <- data.table(
  size = interactions(mats)$size,
  scores = interactions(mats)$score
)
dat <- dat[order(size)]

## Define threshold
thresh <- 300e3

## What perc of loops are >thresh
tb <- table(dat$size > thresh)
(tb[2] / sum(tb)) * 100

## Filter out long (>2Mb) interactions
dat <- dat[size <= thresh]

ns <- round(seq(1, 200, length.out = 9))
par(mfrow=c(3,3))
for (n in ns) {
  dat[, rollMeds := frollapply(scores, n=n, \(x) median(x), align='right')]
  plot(dat$size, dat$rollMeds, type="l", main=paste0("n=", n))
}
par(mfrow=c(1,1))




