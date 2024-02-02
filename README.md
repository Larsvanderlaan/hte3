# hte3 Package

# hte3: Supervised Causal Machine Learning of Heterogeneous Treatment Effects using `sl3`

The `hte3` package, a core component of the Causal Toolbox, equips users with tools for supervised causal machine learning of heterogeneous treatment effects, leveraging the `sl3` package. This package emphasizes precision, robustness, and distribution-free theoretical guarantees.

## Key Features

- **Highly Customizable Meta-learners of HTEs:** Any supervised machine learning algorithm supported by the `sl3` package can be turned into a meta-learner for heterogeneous treatment effects, including the DR-learner, R-learner, T-learner, and EP-learner of the CATE.

- **Novel Meta-learners of the CRATE:** Implements novel EP-learners of the log conditional relative average treatment effect (CRATE), otherwise known as the log conditional relative risk.

- **Data-Adaptive Inverse Weight Truncation:** Automatically perform data-adaptive truncation of estimated inverse propensity weights, enhancing causal inference stability and robustness, especially in scenarios with limited treatment overlap.

- **Generalized Isotonic Calibration:** Calibrate nuisance function estimators, such as outcome regression and propensity score, using generalized isotonic regression. This process enhances the stability, robustness, and bias reduction in causal inferences.

For comprehensive information and practical use cases, consult the package documentation.

### Installation

To install the `hte3` package, use the following command:

```r
devtools::install_github("Larsvanderlaan/causaltoolbox/hte3")
