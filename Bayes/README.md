## Bayesian inference using slicing
In this part of the repo we explore the application and theoretical properties of Slicing in baysian stats.

### Contents: Theory
The goal was to compare performance of MLE and bayesian inference on the logistic regression case on synthetic data.

We simulated data according to the process we want to model:

$$
\begin{align*}
Y_i &\sim Bernoulli(\theta_i) \\
\theta_i &= \text{expit}(\eta_i) \\
\eta_i &= \beta_0+\beta_1 x_i \hspace{1em} x_i\in \mathbb R
\end{align*}
$$

Note how this is only a single dimensional case. Future work could also include more dimensions.

Next we applied MLE and our own function for bayesian logistic regression (`bayes_logit_reg_slice`) and tested it for small samples ranging from 5 to 100. To prove, that the correct prior yields an improvement in performance, we gave the prior the correct beta values. We therefore used the priors:

$$
\begin{pmatrix}
  \tilde \beta_0 \\
  \tilde \beta_1
\end{pmatrix} \sim \mathcal N\left(\begin{pmatrix}
  \beta_0\\
  \beta_1
\end{pmatrix}, \begin{pmatrix}
 2 &0 \\
 0 & 1\\ 
\end{pmatrix}\right)
$$


Finally we tested against 100 test samples with and 20 test runs.

**Results**
- Bayes outperformed by a small margin on very small datasets
- We suspect this gap could increase when we increase the dimension

### Contents: Application
Since we can now apply the function `bayes_logit_reg_slice` to any regression matrix $X$ and dependent variable $y$, we wanted to find a popular dataset and compare performance. To do this, we firstly considered a preprocessed version of the *titanic dataset*, which already was split into train- and test-data.  After application of the logistic regression, we got a test accuracy of 86% which would have put us on place 250 of 16,5k. This enthusiasm sparked the submission we see [here](https://www.kaggle.com/code/carstenstahl/titanic-classification-using-bayes). The final submission only resulted in the 3172th place, which is still pretty good considering there are much more powerfull techniques for binary classification such as XGBoost. 