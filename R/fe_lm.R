#' Apply within transformation and then fit FE regression.
#'
#' @param data A dataframe of genetic data created from sim.phen() or similarly structured
#' @return A regression object.
#' @examples
#' df <- sim.phen()
#' m <- fe.lm(df)

fe.lm <- function(data){
  Y <- data$Y0 - 0.5*(data$Y0 + data$Y1)
  G <- data$g_d0 - 0.5*(data$g_d0 + data$g_d1)
  m <- lm(Y ~ G)
  return(m)
}
