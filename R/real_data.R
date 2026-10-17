## Prepare occupancy (Candanedo) and red-wine (Cortez) copula-scale data
## as in Pan, Nieto-Barajas and Craiu, main.tex Sections 4.4--4.5.

rank01 <- function(x) rank(x, ties.method = "average") / (length(x) + 1)

prepare_occupancy <- function(data_dir) {
  files <- file.path(data_dir, c("datatraining.txt", "datatest.txt", "datatest2.txt"))
  raw <- do.call(rbind, lapply(files, function(f) {
    d <- read.csv(f, check.names = FALSE)
    ## first column may be a row index stored as the name of date
    if (!"date" %in% names(d)) {
      ## datatraining: "date" is column 1 with a leading index column name
      nm <- names(d)
    }
    d
  }))
  ## Harmonise names
  if (!"date" %in% names(raw)) stop("no date column")
  raw$date <- as.POSIXct(raw$date, tz = "UTC")
  day <- as.Date(raw$date, tz = "UTC")
  d5 <- raw[day == as.Date("2015-02-05"), ]
  d5 <- d5[order(d5$date), ]
  ## 5-minute averages of CO2 and HumidityRatio
  d5$slot <- as.POSIXct(trunc(as.numeric(d5$date) / 300) * 300, origin = "1970-01-01", tz = "UTC")
  agg <- aggregate(cbind(CO2, HumidityRatio) ~ slot, data = d5, FUN = mean)
  agg <- agg[order(agg$slot), ]
  ## first differences
  z <- diff(as.matrix(agg[, c("CO2", "HumidityRatio")]))
  U <- apply(z, 2, rank01)
  colnames(U) <- c("CO2", "HR")
  list(U = U, n = nrow(U), tau = cor(U[, 1], U[, 2], method = "kendall"),
       raw_n = nrow(d5), agg_n = nrow(agg))
}

prepare_wine <- function(data_dir) {
  w <- read.csv(file.path(data_dir, "winequality-red.csv"), sep = ";")
  X <- w[, c("fixed.acidity", "citric.acid", "density")]
  U <- apply(as.matrix(X), 2, rank01)
  colnames(U) <- c("fixed_acidity", "citric_acid", "density")
  tau <- cor(U, method = "kendall")
  list(U = U, n = nrow(U), tau = tau)
}
