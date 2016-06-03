\name{dtw}
\alias{dtw}
\title{Dynamic Time Warping}
\usage{
dtw(x,y,window=100)
}
\description{
dtw calculates correspondence of rows of two matrices x and y using dynamic time warping algorithm. Variables x and y should be matrices having the same length of columns. The variable window denotes the calculation window of DTW; if window is large, the result becomes optimum and calculation becomes slower.

Return value is a list.
xsize: number of rows of x.
ysize: number of rows of y.
opt: a two-column matrix that is correspondences of rows of x and y.
}
\examples{

}
