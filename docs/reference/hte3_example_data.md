# Simulate Tiny Example Data for hte3

Convenience generator used in the package examples, tests, and
vignettes.

## Usage

``` r
hte3_example_data(n = 200, seed = 123)
```

## Arguments

- n:

  Number of rows to simulate.

- seed:

  Optional random seed.

## Value

A `data.table` with columns `W1`, `W2`, `W3`, `A`, `Y`, `Y_binary`,
`tau`, `pi1`, `mu0`, and `mu1`.

## Examples

``` r
data <- hte3_example_data(n = 50, seed = 1)
head(data)
#>            W1          W2         W3     A        Y Y_binary        tau
#>         <num>       <num>      <num> <int>    <num>    <int>      <num>
#> 1: -0.4689827 -0.04476076  0.3094479     1 1.735343        1  0.9059141
#> 2: -0.2557522  0.72241895 -0.2936055     0 1.472928        0 -0.7822124
#> 3:  0.1457067 -0.12380579 -0.4594797     0 1.321924        0  0.2245417
#> 4:  0.8164156 -0.51040545  0.9853681     1 5.036641        1  2.9979054
#> 5: -0.5966361 -0.85864191  0.2669865     1 2.593884        0  1.5209631
#> 6:  0.7967794 -0.80106768 -0.5735837     1 4.335661        0  1.4362438
#>          pi1       mu0       mu1
#>        <num>     <num>     <num>
#> 1: 0.4828624 0.9845209 1.8904351
#> 2: 0.3881053 1.4623990 0.6801866
#> 3: 0.4685578 1.5496540 1.7741958
#> 4: 0.7366237 1.9992288 4.9971342
#> 5: 0.5124563 1.2365671 2.7575302
#> 6: 0.5773404 2.4575958 3.8938396
```
