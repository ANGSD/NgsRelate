\name{prune}
\alias{prune}
\title{A function to prune data given a LD tolerance}
\description{
  This function uses the fast and powerfull snpMatrix method of
  calculating LD. 
}
\usage{
prune(data,prune=0.1,LD="D",back=50)
}

\arguments{
  \item{data}{Integer Matrix. A  matrix with SNP genotypes where NA or 0 denotes missing data, 1 for AA, 2 for Aa  and 3 for aa. The number of individuals is the number of rows and the number of SNPs is the number of coloums}
  \item{LD}{The measure used to select the previous SNP to condition on, ("D'","D" or "rsq2")}
  \item{back}{The number of privous SNPs that can be conditioned on (see details for recomadations) }
  \item{prune}{The maximum value allowed for pairwise LD. If 0 then no pruning is performed}
}
\details{

How to select the number of privious SNP that can be conditioned on. First of all if is no LD in the data chosse back=0.If you want to prune SNP away based on LD and/or you want to accomidate LD in the model choose <back> as the number of SNP where you expect there to be LD. For example if you have 500,000 SNP you expect their to be LD between a lot of SNPs in a region. Here I use back=50. If you only have 10,000 SNPs I would use back=5. 
}
\value{
  \item{data}{The pruned dataset}
   \item{usedSnps}{a vector indicating wether a SNP have been used in the analysis (1's) or discarded (0's)}
 
}
\references{ ~put references to the literature/web site here ~ }
\author{Anders Albrechtsen }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
var<-sample(c(NA,1,2,3),10000,replace=TRUE)
data<-matrix(var,ncol=10)
pruned_dataset<-prune(data,0.1,"rsq2",20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
