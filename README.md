# hte3 Package

`This package is in development. Improved documentation and examples coming soon.`
# hte3: Causal Machine Learning of Heterogeneous Treatment Effects using `sl3`

The `hte3` package equips users with tools for supervised causal machine learning of heterogeneous treatment effects, leveraging the `sl3` package. 

## Key Features

- **Highly Customizable Meta-learners of HTEs:** Any supervised machine learning algorithm supported by the `sl3` package can be turned into a meta-learner for heterogeneous treatment effects, including the DR-learner, R-learner, T-learner, and EP-learner of the CATE. For details on the usage of the `sl3` R package, we refer to `https://github.com/tlverse/sl3`.

- **Novel Meta-learners of the CRR:** Implements novel EP-learners of the log conditional relative risk (CRR).


For comprehensive information, consult the package documentation.

### Installation

To install the `hte3` package, use the following command:

```r
if(!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("Larsvanderlaan/hte3")
```
