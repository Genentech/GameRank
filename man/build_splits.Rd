% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_splits.R
\name{build_splits}
\alias{build_splits}
\title{Helper function to generate splits that are likely to produce models}
\usage{
build_splits(m, dat, resp, vars, fn_train, fn_eval, min_cc = 25L, ...)
}
\arguments{
\item{m}{Number of parallel random splits to generate}

\item{dat}{data.frame or tibble comprising data for model generation and validation.}

\item{resp}{Character string defining the response (lhs) of the model formula.}

\item{vars}{Character vector defining the list of variables for selection. Those are concatenated by '+' 
as the right hand side (rhs) of the modelling formula.}

\item{fn_train}{Function with signature function( dat, resp, selection, ... ) that returns a model or NULL in any other case on the given data dat.}

\item{fn_eval}{Function with signature function( dat, resp, selection, ... ) that returns a real number or NA in any other case, e.g. when model is NULL.}

\item{min_cc}{Minimal number of complete case observations to continue with building splits}

\item{...}{An other arguments passed to fn_train or fn_eval during calls, e.g. maybe 'u = 365' for Survival evaluations specifying the landmark day.}
}
\value{
Matrix with nrow(dat) rows and m columns comprising 1s and 2s to indicate if sample is part of training or validation split.
}
\description{
The function should generate m random splits for dat, resp, vars such
that fn_train is likely to produce models. E.g. if there are a lot of missing values
in the dataset and splitting reduces these in an unfavorable way, model parameters for lm
or glm might not be estimable. That happens likely if the number of complete.cases is too low.
Therefore this function first separates the complete and incomplete cases, and then
randomly samples and distributes them into training and validation splits. This ensures that
both splits include 2/3rds and 1/3rd complete cases each.
}
\examples{
resp <- "resp"
vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )

ds <- build_splits( 3L, toy_data, resp, vars, fn_train_binomial, fn_eval_binomial_auroc )
ds

}
