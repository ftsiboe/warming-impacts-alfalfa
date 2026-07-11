#' Apply the delta method to a list of formulas or expressions
#'
#' @description
#' Evaluates a list of user-supplied formulas or expressions using
#' `car::deltaMethod()` and returns a tidy data frame of transformed estimates
#' and associated inference results.
#'
#' This is useful when you want to compute nonlinear functions of estimated
#' coefficients and obtain standard errors and confidence intervals using the
#' delta method.
#'
#' @param func A list of formulas or expressions defining the transformations to
#'   evaluate. Each element should contain a right-hand-side expression that can
#'   be interpreted by `car::deltaMethod()`.
#' @param vcMat A variance-covariance matrix for the coefficient estimates.
#'   Defaults to `VCOV`.
#' @param coefs A named numeric vector of coefficient estimates. Defaults to
#'   `COEF`.
#'
#' @return A data frame with one row per transformation. The returned object
#'   includes:
#' \describe{
#'   \item{form}{Character representation of the transformation.}
#'   \item{Estimate and inferential columns}{Output returned by
#'   `car::deltaMethod()`, excluding the last two columns produced internally.}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Extracts the right-hand-side expression from each element of `func`.
#'   \item Converts each expression to character form.
#'   \item Applies `car::deltaMethod()` to each transformation using the supplied
#'   coefficients and variance-covariance matrix.
#'   \item Combines the results into a single data frame.
#'   \item Removes the final two columns of the intermediate result before
#'   returning the output.
#' }
#' @export
compute_delta_method <- function(func, vcMat = VCOV, coefs = COEF) {
  # func  <- Form
  # coefs <- coef(PWM)
  # vcMat <- PWM$vcov

  temp <- lapply(func, function(x) as.character(as.expression(x[[length(x)]])))
  func <- data.frame(form = unlist(temp))

  lisRes <- apply(
    func,
    1,
    function(x) car::deltaMethod(
      object = coefs,
      g = x,
      vcov. = vcMat,
      level = 0.99
    )
  )

  val <- plyr::ldply(lisRes)
  val <- cbind(func, val)
  lenVal <- length(val[1, ])
  retDF <- val[, c(-(lenVal - 1), -lenVal)]

  return(retDF)
}
