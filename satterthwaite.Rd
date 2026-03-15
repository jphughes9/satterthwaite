\name{satterthwaite.lmerMod, satterthwaite.glmerMod}
\alias{satterthwaite}

\title{
Compute Satterthwaite degrees of freedom for a hypothesis test for an lmerMod/glmerMod object
}
\description{
\code{satterthwaite} returns denominator degrees of freedom for a test
of the hypothesis \deqn{ Ho: L \beta = 0} for lmer/glmer objects, where \eqn{L} is 
a \eqn{q \times p} matrix of contrasts and \eqn{\beta} is the vector of fixed effects.
}
\usage{
satterthwaite(obj,L,simplify=FALSE)
}

\arguments{
 \item{obj}{
Fitted lmer/glmer object
}
\item{L}{
Matrix of (independent) contrasts of dimension q x p where q is the number contrasts and p is the 
number of fixed effects
}
\item{simplify}{
If one or more (but not all) of the variance components are zero on the theta scale,
\code{simplify=TRUE} tells satterthwaite to compute the degrees of freedom after 
eliminating the zero variance components from the calculation (and keeping the 
same fixed effects). If \code{simplify=FALSE} then the function returns an error if there 
are any zero variance components (on the theta scale). Default is \code{simplify=FALSE}.
Note that the function returns an error when all variance components are zero regardless
of the value of \code{simplify}.
}
}

\details{
The lmerMod function implements the df_1d and df_md function from the mmrm package for lmer models.
See accompanying article (Satterthwaite’s Method for Degrees of Freedom in Linear Mixed Models by Rune Haubo B Christensen)
for details on the computational methods. The glmerMod function makes some changes, primarily by computing the 
variance-covariance matrix of the parameters using the internally available Hessian (\code{model@optinfo$derivs$Hessian}) 
 and derivatives of the fixed effect covariance matrix
the internal function \code{pp$unsc()} (rather than the \code{pp$RXi()} function).

In satterthwaite.glmerMod(), if there is a variance component on the boundary
(defined as <.001 in absolute value on the theta scale), then a warning is issued
 and the degrees of freemdom calculation is done on (the
equivalent of) a simpler model without that variance component. If all variance 
components are on the boundary, the function returns an error.

The satterthwaite function gives the 
denominator degrees of freedom for an F-test of the hypothesis Ho. A p-value is computed as
\code{1-pf(f,q,df)} where \code{f} is the computed f-statistic. For a single contrast the
p-value is computed as
\code{1 - pt(t,df)} where \code{t} is the computed t-statistic.
}

\value{
For 1 degree of freedom tests, a list with elements:

\item{estimate}{Point estimate of the contrast}
\item{se}{Standard error of the contrast}
\item{tstat}{t-statistic for the test Ho}
\item{ddf}{Satterthwaite degrees of freedom for the t-test of Ho}
\item{pvalue}{p-value for the t-test}


For multi-degree of freedom tests, a list with elements:
\item{F value}{F-statistic for the test Ho}
\item{q}{Numerator degrees of freedom for the F-test of Ho}
\item{ddf}{Satterthwaite denominator degrees of freedom for the F-test of Ho}
\item{pvalue}{p-value for the F-test}
}

\references{
Bates, D et al. Fitting Linear Mixed-Effects Models Using lme4. J of Staistical Software, 67(1), 2015. doi: 10.18637/jss.v067.i01

Christensen, RHB. Satterthwaite’s Method for Degrees of Freedom in Linear Mixed Models. Posted with the mmrm R package

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
L = matrix(c(0,0,0,0,1,0,-1,0,0,0,0,1,-1,0),2,7,byrow=TRUE)
satterthwaite(rslt,L)
}
\keyword{ satterthwaite degrees of freedom }
\keyword{ mixed models }% __ONLY ONE__ keyword per line
