#' Is InChIKey
#'
#' Check if a string is a validly-formatted InChIKey
#'
#' @details The criteria are the following, as described in Heller et al. (2015)
#'   (see References).
#'
#' - must have 27 characters
#' - the first 14 characters must be uppercase letters
#' - the 15th character must be a hyphen
#' - the 16th-23rd characters must be uppercase letters
#' - the 24th character must be either "S" or "N"
#' - the 25th character must be "A"
#' - the 26th character must be a hyphen
#' - the 27th character must be an uppercase letter
#'
#'
#' @param x Character: A vector of one or more strings to check for InChIKey
#'   format.
#' @return A logical vector the same length as `x`, containing `TRUE` if valid
#'   InChIKey format, and `FALSE` otherwise.
#' @export
#' @author Caroline Ring
#' @references Heller, S.R., McNaught, A., Pletnev, I. et al. InChI, the IUPAC
#'   International Chemical Identifier. J Cheminform 7, 23 (2015).
#'   https://doi.org/10.1186/s13321-015-0068-4

is_inchikey <- function(x){
  if(is.character(x)){
    return(
      grepl(x=x,
          pattern = "^[[:upper:]]{14}\\-[[:upper:]]{8}[SN]A\\-[[:upper:]]$")
    )
    #must have 27 characters
  #the first 14 must be uppercase letters only
    #the 15th must be a hyphen
    #the 16-23rd must be uppercase letters only
    #the 24th must be either S or N
    #the 25th must be A
    #the 26th must be a hyphen
    #the 27th must be an uppercase letter
  }else{
    return(rep(FALSE, length(x)))
  }
}
