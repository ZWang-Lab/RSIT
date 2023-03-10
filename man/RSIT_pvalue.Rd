% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RSIT_pvalue.R
\name{RSIT_pvalue}
\alias{RSIT_pvalue}
\title{Calculate  P-values}
\usage{
RSIT_pvalue(TYall, G, impute.method = "fixed")
}
\arguments{
\item{TYall}{The output of function "RSIT_est()"}

\item{G}{The genotype matrix, an m*q matrix where m is the number of subjects and q is the total number genetic variants.}

\item{impute.method}{choose the imputation method when there is missing genotype. Optional options are: 'random', 'fixed' or 'bestguess'.}
}
\value{
This function returns a vector of the P-value of RSIT-B, RSIT-S, and RSIT-E.
}
\description{
This function tests a SNPs for a given SNP set for a given RSIT estimated null model.
}
