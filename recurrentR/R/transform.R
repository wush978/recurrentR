#'@title Transformate Data to the \code{recurrent-data} Object
#'@param src data.frame. This is the raw data frame.
#'@param id character value. The column name in \code{src} for \emph{id}.
#'@param time character value.
#'@param time_type character value.
#'@param indicator character value.
#'@param indicator_value list of named value.
#'@param covariate character value.
#'@return S4 Object of \code{recurrent-data}
#'@export
#'@examples
#'\dontrun{
#'library(survrec)
#'data(MMC)
#'obj <- create_recurrent_data.data.frame(
#'  MMC, id = "id", time = "time", time_type = "relatively",
#'  indicator = "event", indicator_value = list("recurrent" = 1, "censor" = 0),
#'  covariate = "group"
#')
#'}
create_recurrent_data.data.frame <- function(src, id, time, time_type = c("absolutely", "relatively"),
                                             indicator, indicator_value, covariate) {
  spec <- list(id = id, time = time, time_type = time_type, indicator = indicator, indicator_value = indicator_value, covariate = covariate)
  id.index <- split(seq_len(nrow(src)), src[[spec$id]])
  id.group <- lapply(id.index, function(i) src[i,])
  y <- as.vector(sapply(USE.NAMES=FALSE, id.group, function(df) {
    stopifnot(df[-nrow(df),spec$indicator] == spec$indicator_value$recurrent)
    stopifnot(df[nrow(df),spec$indicator] == spec$indicator_value$censor)
    switch(
      spec$time_type,
      "absolutely" = {
        stopifnot(df[[spec$time]] >= 0)
        stopifnot(diff(df[[spec$time]]) >= 0)
        df[nrow(df),spec$time]
      },
      "relatively" = {
        stopifnot(df[[spec$time]] >= 0)
        sum(df[[spec$time]])
      },
      stop("Invalid time_type")
    )
  }))
  T_0 <- if(is.null(spec$T_0)) max(y) else spec$T_0
  D <- y < T_0
  t <- lapply(id.group, function(df) {
    switch(
      spec$time_type,
      "absolutely" = {
        head(df[[spec$time]], -1)
      },
      "relatively" = {
        head(cumsum(df[[spec$time]]), -1)
      }
    )
  })
  W <- do.call(rbind, lapply(id.group, function(df) {
    df[1,spec$covariate]
  }))
  create_recurrent_data.numeric(y, D, t, T_0, W)
}
