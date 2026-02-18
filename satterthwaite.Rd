\name{satterthwaite.lmerMod}
\alias{satterthwaite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Sattertheaite degrees of freedom for a hypothesis test for an lmerMod object
}
\description{
\code{satterthwaite} returns denominator degrees of freedom for a test
of the hypothesis \deqn{ Ho: L \beta = 0} for lmer objects, where \eqn{L} is 
a \eqn{r \times p} matrix of contrasts and \eqn{\beta} is the vector of fixed effects.
}
\usage{
satterthwaite(obj,lvec)
}

\arguments{
 \item{obj}{
Fitted lmer object
}
\item{lvec}{
Matrix of contrasts of dimension r x p where r is the number contrasts and p is the 
number of fixed effects
}
}

\details{
See accompanying article for details on the computational methods. The function gives the 
denominator degrees of freedom (possibly 0) for an F-test of the hypothesis Ho. A p-value may be computed as
\code{1-pf(f,r,df)} where \code{f} is the computed f-statistic. For a single contrast one may equivalently compute
\code{1 - pt(t,df)} where \code{t} is the computed t-statistic.
}

\value{
A real number (possibly 0) giving the denominator degrees of freedom for an 
F-test of the hypothesis Ho.
}

\references{
Bates, D et al. Fitting Linear Mixed-Effects Models Using lme4. J of Staistical Software, 67(1), 2015. doi: 10.18637/jss.v067.i01

SAS GLIMMIX documentation (https://documentation.sas.com/doc/en/statug/15.2/statug_glimmix_syntax01.htm)
}
\author{
James P Hughes
}
%%\note{
%%  ~~further notes~~
%%}

\examples{
rslt = lmer(outvar ~ ftime + ftreat + (1 | fcluster), data=data)
lvec = matrix(c(0,0,0,0,1,0,-1,0,0,0,0,1,-1,0),2,7,byrow=TRUE)
satterthwaite(rslt,lvec)
}
\keyword{ satterthwaite degrees of freedom }
\keyword{ mixed models }% __ONLY ONE__ keyword per line
