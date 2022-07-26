# GameRank - An R package for Feature Selection

The GameRank package is an R package that implements functions for 
variable and feature selection to create predictive models. This includes a 
number of feature selection algorithms, mostly combinatorial
search strategies, which are:

 * Forward selection,
 * Backward selection,
 * Bidirectional search, and
 * Plus-L, Minus-R search together with 
 * Random search, which draws random combinations across the feature space, and
 * *GameRank* a novel feature selection algorithm, that employs a
   maximum likelihood ranking model to learn the optimal ranking of features. 
   The following figure illustrates the algorithm: 
   ![The GameRank Algorithm](man/figures/GameRank-Algorithm.png)
   
In addition, it provides methods for feature construction, that is

 * functions to screen variable properties,
 * evaluate simple variable transformations, 
 * Box-Cox transformations for regression and binomial models, as well as
 * methods for the detection of features with multi-modal distributions.
 
GameRank allows for direct optimization of measures related to _calibration_ and
_discrimination_ of models. This approach uses a training:validation:test split-set
approach. However, for small sample selection projects, functions for bias-corrections via
cross-validation or bootstrapping are included.

Parts of this code have been used to craft a model for chemo-tolerability 
prediction. [https://ascopubs.org/doi/full/10.1200/CCI.21.00121]

The package can be installed directly from GitHub via

```{r}
devtools::install_github("Genentech/GameRank")
```

Here are some few lines of example code to use the package to 
**screen features:**
```{r}
library( dplyr )
library( ggplot2 )
vck <- check_variables( toy_data, response_var, list_variables )
vck %>% 
   ggplot(aes(x=entropy, y=mutual_information) ) +
   geom_point()
```
and do a **forward selection:**
```{r}
sel <- forward( toy_data, response_var, list_variables 
                 fn_train_binomial,  # Train a logistic regression model
                 fn_eval_binomial_auroc,  # Evaluate the Area under ROC curve
                 m=4L,   # Search for combination of 4 features
                 splits=2L, # Evaluate on 2 parallel split sets
                 maximize=TRUE # Maximize the Area under ROC curve
                 )
sel$variable_selections
```

Have fun!