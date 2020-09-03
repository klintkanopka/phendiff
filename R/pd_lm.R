#' Fit PD regression.
#'
#' @param data A dataframe of genetic data created from sim.phen() or similarly structured
#' @param rho Genetic correlation between siblings. Defaults to 0.5
#' @return A regression object.
#' @examples
#' df <- sim.phen()
#' m <- pd.lm(df, 0.5)

pd.lm <- function(data, rho=0.5){
  dY <- data$Y0 - data$Y1
  G <- (1-rho) * data$g_d0
  m <- lm(dY ~ G)
  return(m)
}
