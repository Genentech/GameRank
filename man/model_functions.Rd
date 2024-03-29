% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_functions.R, R/model_functions_ml_imp.R
\name{model_functions}
\alias{model_functions}
\alias{fn_train_normal}
\alias{fn_eval_normal}
\alias{fn_train_binomial}
\alias{fn_eval_binomial_mad}
\alias{fn_eval_binomial_auroc}
\alias{fn_train_cox}
\alias{fn_eval_cox}
\alias{fn_train_lda}
\alias{fn_eval_lda}
\alias{fn_train_qda}
\alias{fn_eval_qda}
\alias{fn_train_normal_ml_imp}
\alias{fn_train_binomial_ml_imp}
\alias{fn_train_cox_ml_imp}
\title{Functions for Model fitting and evaluation}
\usage{
fn_train_normal(dat, resp, selection, ...)

fn_eval_normal(dat, resp, selection, mod, ...)

fn_train_binomial(dat, resp, selection, ...)

fn_eval_binomial_mad(dat, resp, selection, mod, ...)

fn_eval_binomial_auroc(dat, resp, selection, mod, ...)

fn_train_cox(dat, resp, selection, ...)

fn_eval_cox(dat, resp, selection, mod, u = NULL, ...)

fn_train_lda(dat, resp, selection, lda_fit_type = "mle", ...)

fn_eval_lda(dat, resp, selection, mod, lda_pred_type = "plug-in", ...)

fn_train_qda(dat, resp, selection, qda_fit_type = "mle", ...)

fn_eval_qda(dat, resp, selection, mod, qda_pred_type = "plug-in", ...)

fn_train_normal_ml_imp(dat, resp, selection, ...)

fn_train_binomial_ml_imp(dat, resp, selection, ...)

fn_train_cox_ml_imp(dat, resp, selection, ...)
}
\arguments{
\item{dat}{data.frame or tibble rows from the full dataseet provided to the
wrapper that should be used for generating or evaluating models.}

\item{resp}{Response variable being the lhs of the model formula}

\item{selection}{Current selection for model generation or evaluation}

\item{...}{Any other arguments passed to both types of functions, 
e.g. 'u = 365' to define the landmark day for survival probability evaluation.}

\item{mod}{For evaluation functions the model to be evaluated on dat}

\item{u}{Landmark time point at which the survival probability is evaluated}

\item{lda_fit_type}{Type parameter for MASS::lda function, which model to fit}

\item{lda_pred_type}{Prediction method parameter for MASS::predict.lda}

\item{qda_fit_type}{Type parameter for MASS::qda function, which model to fit}

\item{qda_pred_type}{Prediction method parameter for MASS::predict.qda}
}
\value{
fn_train_... functions return a fitted model object or NULL if that 
fails. fn_eval_... functions return a numeric real value measuring the 
validation performance on the given data, or NA if that fails.
}
\description{
The interface for model generating functions is
\strong{function( dat, resp, selection, ... )}
where dat is a data.frame or tipple sliced out of the overall
dataset according to the split matrix (1 = training) and (2 = validation).
Furthermore, resp denotes a character string that should be used as
response variable (on the lhs of a formula) and selection denotes
a vector of character strings that represents the feature combination to
be evaluated.
Model function should return NULL if the model cannot be fit or determined or
other reasons exist for why it cannot be generated.

The interface for evaluation functions is \strong{function( dat, resp, 
selection, mod, ... )} with the same parameters as model generating functions 
but extended by the model mod. Evaluation functions should return a numeric 
value and NA in case of the model is NULL or any other failure preventing
evaluation.

Model functions ending with _imp modify the selection variables to perform a
maximum-likelihood imputation for missing values. Precisely, let \eqn{x} be 
the variable it is replaced by \eqn{\phi(x) + \psi(x)} where 
\eqn{\phi(x) = x if x is not missing, and 0 otherwise} and
\eqn{\psi(x) = 1 if x is missing, and 0 otherwise}. This equates to 
a 0-1-Dummy coding for missing values and the effect for \eqn{\psi} can be
interpreted as the maximum-likelihood imputation value to be used if \eqn{x} 
is missing.
}
