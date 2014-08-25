#-----------------------------------------------------------------------
#     Copyright (C) 2004-2013  Serge Iovleff
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
#' rtkpp is a R/stk++ bridge package based on the existing STK++ library.
#'
#' This package contains the header files for the STK++ template
#' library. The typical usage is to install this package and list it in
#' the \env{LinkingTo: } line in the \file{DESCRIPTION} file of
#' other packages.
#'
#' As described at the STK++ project's home page, \url{http://www.stkpp.org},
#' STK++ is a versatile, fast, reliable and elegant collection of C++ classes
#' for statistics, clustering, linear algebra, arrays (with  Eigen),
#' regression, dimension reduction, etc. Some functionalities provided by the
#' library are available in the R environment.
#'
#' \tabular{ll}{
#'   Package: \tab rtkpp\cr
#'   Type: \tab Package\cr
#'   Version: \tab 2.0.1\cr
#'   Date: \tab 2014-07-05\cr
#'   License: \tab GPL for the rtkpp side, LGPL for the stkpp side  + file LICENSE\cr
#'   LazyLoad: \tab yes\cr
#' }
#'
#' @rdname rtkpp-package
#' @name rtkpp
#' @aliases rtkpp
#' @docType package
#' @keywords stk++, stkpp
#' @import Rcpp
#'
#' @author
#' Author: Serge Iovleff \email{contact@@stkpp.org}
#'
#' @useDynLib rtkpp
NULL
#' Quantitative data: Old Faithful Geyser
#'
#' The file geyser.rda contains 272 observations from the Old Faithful Geyser
#' in the Yellowstone National Park. Each observation consists of two
#' measurements: the duration (in minutes) of the eruption and the waiting
#' time (in minutes) to the next eruption.
#'
#' Old Faithful erupts more frequently than any other big geyser, although it
#' is not the largest nor the most regular geyser in the park. Its average
#' interval between two eruptions is about 76 minutes, varying from
#' 45 - 110 minutes. An eruption lasts from 1.1/2 to 5 minutes,
#' expels 3,700 - 8,400 gallons (14,000 - 32,000 liters) of boiling water, and
#' reaches heights of 106 - 184 feet (30 - 55m). It was named for its consistent
#' performance by members of the Washburn Expedition in 1870. Old Faithful is
#' still as spectacular and predictable as it was a century ago.
#'
#' @format A data frame with 272 observations on the following 2 variables.
#'
#' \describe{
#'   \item{\code{Duration}}{a numeric vector containing the duration (in minutes) of the eruption}
#'   \item{\code{Waiting.Time}}{a numeric vector containing the waiting time (in minutes) to the next eruption}
#' }
#'
#' @source \url{http://www.geyserstudy.org/geyser.aspx?pGeyserNo=OLDFAITHFUL}
#'
#' @references
#' Hardle, W. (1991). "Smoothing Techniques with Implementation in S". Springer-Verlag, New York.
#' Azzalini, A. and Bowman, A. W. (1990). "A look at some data on the Old Faithful geyser". Applied Statistics 39, 357-365.
#'
#' @name geyser
#' @docType data
#' @keywords datasets
#'
#' @examples
#'   data(geyser)
NULL
