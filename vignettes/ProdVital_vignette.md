---
title: "Estimating growth increment, productivity, and vital rates with ProdVital"
output: 
  html_document:
    toc: true
author: "John L. Godlee"
date: 2023-08-10
bibliography: ProdVital_vignette.bib
---



# Setup

Load the `{seosawr}` package, which contains some useful functions for cleaning
SEOSAW stem data. The `{seosawr}` package is located in the [SEOSAW Bitbucket
repository](https://bitbucket.org/miombo/seosaw). Refer to the `README.md` in
that repository for how to install the package. Load the `{dplyr}` package for
basic data manipulation. Load `{ProdVital}` which contains the core functions
for analysing phyto-demographic data.


```r
# Packages
library(dplyr)
library(seosawr)
library(ProdVital)
library(ggplot2)
```

Load the example dataset included in `{ProdVital}`:


```r
data(bicuar)
```

`bicuar` is a dataframe of stem data from two plots in Bicuar National Park,
southwest Angola, measured across three censuses. Censuses were conducted in
2019, 2021 and 2023. The original data is real, but it has been adjusted to
better illustrate the functionality of `{ProdVital}`. The data is in the SEOSAW
style, but with columns removed where they are not relevant to this workflow.
Each row is a stem measurement in a given census.

* `plot_id` - ID of the plot
* `census_date` - Census date of measurement 
* `stem_id` - ID value grouping measurements into stems
* `stem_status` - Vital status: alive, resprouting, dead
* `diam` - Stem diameter (cm)
* `agb` - Estimated aboveground biomass (Mg)

# Prepare the stem data

Stems may be missing from a particular census, despite being measured in
censuses before and after. Stem timelines should not have gaps, otherwise you
will get unrealistic estimates of recruitment, mortality, and biomass turnover
using the functions in `{ProdVital}`. We need to create these rows so we can
fill them with estimated values. We can do this using the
`seosawr::stemFillGen()` function:


```r
# Generate missing rows
plots <- unique(bicuar[, c("plot_id", "census_date")])
missing_records <- stemFillGen(bicuar, plots)

# All the new records should be in the middle census
stopifnot(all(missing_records$census_date == 2021))

# Add the missing rows to the stem data
bicuar_pad <- bind_rows(bicuar, missing_records)
```

Sometimes, a stem will appear to have died only to resurrect at a later date.
Stems might also be missing stem status values. We can try to fix these issues
using the `seosawr::statusImputGen()` function. 


```r
# Impute mortality status
status_est <- statusImputGen(bicuar_pad)

# Add new columns to stem data
bicuar_sfix <- bind_cols(bicuar_pad, status_est)
```

Now that the stem status values have been fixed, we should remove measurements
of dead stems. The functions in `{ProdVital}` assume that the stem data only
contain measurements from living stems. Measurements from dead stems should be
excluded from the stem data. The absence of a stem from a particular census
after it was initially measured implies that the stem was dead during that
census. Similarly the absence of a stem from a particular census before it was
initially measured implies that the stem had not yet recruited. It's up to you
whether you count resprouting stems (`stem_status` "r", stem dead above point
of measurement but with living material elsewhere on the tree) as alive or
dead:


```r
bicuar_alive <- bicuar_sfix %>%
    filter(stem_status_est == "a")
```

Stem diameter measurements might be missing. These values should be imputed if
you are calculating plot-level productivity, otherwise you will get unrealistic
estimates of growth and loss. A simple way to impute missing diameter values is
to take the mean of measured diameters before and after the missing diameter:


```r
bicuar_clean <- bicuar_alive %>%
    group_by(plot_id, stem_id) %>%
    mutate(diam = ifelse(is.na(diam), mean(diam, na.rm = TRUE),
        diam), agb = ifelse(is.na(agb), mean(agb, na.rm = TRUE),
        agb))

stopifnot(all(!is.na(bicuar_clean$diam)))
```

It is important to understand that this naive approach doesn't consider that
some diameter measurements may contain errors, what should happen if the
missing value is at the start or end of the stem timeline, or what should
happen if there are multiple consecutive missing values. In these cases, more
complex approaches will be required to impute missing diameter values.

If you only want to calculate stem growth increment, these missing values can
simply be removed:


```r
bicuar_excl <- bicuar_clean %>%
    filter(!is.na(diam))
```

# Identify stem fate

Using the `obsID()` function, it is possible to identify the fate of individual
stems within a given census interval. __Recruits__ are stems which were not
present in the initial census, but are present in the final census.
__Survivors__ are stems which are present in both the initial and final census.
__Deaths__ are stems which are present in the initial census but not in the
final census:


```r
bicuar_clean$plot_stem <- paste(bicuar_clean$plot_id, bicuar_clean$stem_id,
    sep = ":")

obs1 <- obsID(bicuar_clean, t0 = "2019", tT = "2021", type = c("rec",
    "sur", "mor"), group = "plot_stem", census = "census_date")

lapply(obs1, head)
```

```
## $rec
##      plot_stem
## 2067 ABG_5:679
## 2075 ABG_5:680
## 2077 ABG_5:681
## 2079 ABG_5:682
## 2080 ABG_5:683
## 2082 ABG_5:684
## 
## $sur
##     plot_stem
## 1     ABG_5:1
## 4    ABG_5:10
## 7   ABG_5:100
## 10 ABG_15:101
## 12  ABG_5:101
## 16  ABG_5:102
## 
## $mor
##      plot_stem
## 83   ABG_5:114
## 118 ABG_15:121
## 145 ABG_15:126
## 164   ABG_5:13
## 520  ABG_5:203
## 587  ABG_5:226
```

Note that `plot_stem` was created to provide single ID column, but this is also
valid:


```r
obs2 <- obsID(bicuar_clean, t0 = "2019", tT = "2021", type = c("rec",
    "sur", "mor"), group = c("plot_id", "stem_id"), census = "census_date")

lapply(obs2, head)
```

```
## $rec
##      plot_id stem_id
## 2067   ABG_5     679
## 2075   ABG_5     680
## 2077   ABG_5     681
## 2079   ABG_5     682
## 2080   ABG_5     683
## 2082   ABG_5     684
## 
## $sur
##    plot_id stem_id
## 1    ABG_5       1
## 4    ABG_5      10
## 7    ABG_5     100
## 10  ABG_15     101
## 12   ABG_5     101
## 16   ABG_5     102
## 
## $mor
##     plot_id stem_id
## 83    ABG_5     114
## 118  ABG_15     121
## 145  ABG_15     126
## 164   ABG_5      13
## 520   ABG_5     203
## 587   ABG_5     226
```

`obsID()` serves as a building block in many other functions in `{ProdVital}`,
like `prodKohyama()`, `vitalKohyama()` and `obsSurGrowth()`.

# Observed stem growth, recruitment growth, and mortality loss

The `obsSurGrowth()` function can calculate the growth increment of all stems
which survived between two censuses. The argument `w` gives the column used to
calculate growth. The function returns a named vector of stem growth increments
between the two census dates. The names refer to the `group` ID column(s) used
to define individuals. 


```r
obsSurGrowth(bicuar_clean, t0 = "2019", tT = "2021", w = "diam",
    group = c("plot_id", "stem_id"), census = "census_date")
```

```
##    ABG_5.1   ABG_5.10  ABG_5.100 ABG_15.101  ABG_5.101  ABG_5.102 ABG_15.103 
##      0.400      1.000      0.000      0.200     -0.100     -0.100     -0.200 
##  ABG_5.103 ABG_15.104  ABG_5.104 ABG_15.105  ABG_5.105 ABG_15.106  ABG_5.106 
##     -0.800     -0.700     -0.400      0.000     -0.800      0.200     -0.500 
## ABG_15.107  ABG_5.107 ABG_15.108  ABG_5.108 ABG_15.109  ABG_5.109   ABG_5.11 
##      0.500     -0.500      0.000     -0.100     10.700     -2.900      0.200 
## ABG_15.110  ABG_5.110 ABG_15.111  ABG_5.111 ABG_15.112  ABG_5.112 ABG_15.113 
##     -0.500     -0.500      0.400      0.600      0.400      1.300      0.500 
##  ABG_5.113 ABG_15.114 ABG_15.115  ABG_5.115 ABG_15.116  ABG_5.116 ABG_15.117 
##      0.000      0.000     -0.100     -0.300      0.000     -0.100      0.000 
##  ABG_5.117 ABG_15.118  ABG_5.118 ABG_15.119  ABG_5.119   ABG_5.12 ABG_15.120 
##      0.200     -0.200     -0.400      0.200      2.900     -0.200      4.200 
##  ABG_5.120  ABG_5.121 ABG_15.122  ABG_5.122 ABG_15.123  ABG_5.123 ABG_15.124 
##     -0.200      0.000      0.800      0.700     -0.100     -1.300     -0.800 
##  ABG_5.124 ABG_15.125  ABG_5.125  ABG_5.126 ABG_15.127  ABG_5.127 ABG_15.128 
##     -0.900      0.000     -0.200     -0.600     -0.600      1.300      0.000 
##  ABG_5.128 ABG_15.129  ABG_5.129 ABG_15.130  ABG_5.130 ABG_15.131  ABG_5.131 
##     -0.300      2.043      0.200      0.500      0.600      0.000      0.100 
## ABG_15.132  ABG_5.132  ABG_5.133 ABG_15.134  ABG_5.134 ABG_15.135  ABG_5.135 
##      0.100     -0.500      0.300     -1.200     -0.100     -0.200     -1.000 
## ABG_15.136  ABG_5.136  ABG_5.137 ABG_15.138  ABG_5.138 ABG_15.139  ABG_5.139 
##      0.000     -0.400     -0.300      0.000      0.100      0.500      2.361 
##   ABG_5.14 ABG_15.140  ABG_5.140 ABG_15.141  ABG_5.141 ABG_15.142  ABG_5.142 
##      0.200     -0.100      0.800      0.000      0.000      0.000      0.400 
## ABG_15.143  ABG_5.143 ABG_15.144  ABG_5.144 ABG_15.145  ABG_5.145 ABG_15.146 
##      0.000      0.500      0.000      0.000     -0.200     -0.100      0.400 
##  ABG_5.146 ABG_15.147  ABG_5.147 ABG_15.148  ABG_5.148 ABG_15.149  ABG_5.149 
##      0.000      0.000      1.000      0.500      1.000      0.100      4.000 
##   ABG_5.15 ABG_15.150 
##      0.100      0.300 
##  [ reached getOption("max.print") -- omitted 693 entries ]
```

With options `full = TRUE`, `obsSurGrowth()` returns a dataframe rather than a
named vector, which may be easier to join back to the original data using the
columns in `group`.


```r
obsSurGrowth(bicuar_clean, t0 = "2019", tT = "2021", w = "diam",
    group = c("plot_id", "stem_id"), census = "census_date",
    full = TRUE)
```

```
##     plot_id stem_id obs_sur_growth
## 299   ABG_5       1            0.4
## 302   ABG_5      10            1.0
## 303   ABG_5     100            0.0
## 2    ABG_15     101            0.2
## 306   ABG_5     101           -0.1
## 308   ABG_5     102           -0.1
## 4    ABG_15     103           -0.2
## 309   ABG_5     103           -0.8
## 6    ABG_15     104           -0.7
## 311   ABG_5     104           -0.4
## 7    ABG_15     105            0.0
## 314   ABG_5     105           -0.8
## 10   ABG_15     106            0.2
## 316   ABG_5     106           -0.5
## 11   ABG_15     107            0.5
## 317   ABG_5     107           -0.5
## 13   ABG_15     108            0.0
## 320   ABG_5     108           -0.1
## 16   ABG_15     109           10.7
## 322   ABG_5     109           -2.9
## 323   ABG_5      11            0.2
## 17   ABG_15     110           -0.5
## 325   ABG_5     110           -0.5
## 19   ABG_15     111            0.4
## 327   ABG_5     111            0.6
## 22   ABG_15     112            0.4
## 329   ABG_5     112            1.3
## 23   ABG_15     113            0.5
## 331   ABG_5     113            0.0
## 26   ABG_15     114            0.0
## 27   ABG_15     115           -0.1
## 333   ABG_5     115           -0.3
## 30   ABG_15     116            0.0
##  [ reached 'max' / getOption("max.print") -- omitted 760 rows ]
```

Similarly, the `obsMorLoss()` function calculates the loss from stems which
died within a census interval. The arguments of the function are identical to
`obsSurGrowth()`:


```r
obsMorLoss(bicuar_clean, t0 = "2019", tT = "2021", w = "diam",
    group = c("plot_id", "stem_id"), census = "census_date")
```

```
## ABG_15.121 ABG_15.126  ABG_5.114   ABG_5.13  ABG_5.203  ABG_5.226  ABG_5.228 
##       28.5        5.0       55.4        5.0       58.5        6.5        5.5 
##  ABG_5.239   ABG_5.24  ABG_5.256  ABG_5.261  ABG_5.322  ABG_5.326  ABG_5.328 
##        5.2        6.1        9.0        9.2        8.0        5.2       12.4 
##  ABG_5.335  ABG_5.339  ABG_5.343  ABG_5.345  ABG_5.347  ABG_5.360  ABG_5.361 
##        7.0        7.9        5.2        6.0        9.9       10.3        6.6 
##  ABG_5.362  ABG_5.363  ABG_5.366  ABG_5.368  ABG_5.392  ABG_5.393  ABG_5.406 
##        7.5        6.5        9.5        6.0        7.8        7.7       10.0 
##  ABG_5.422  ABG_5.453  ABG_5.464  ABG_5.475   ABG_5.49  ABG_5.540  ABG_5.549 
##       11.7        8.2        7.9        6.0        5.0        8.3        5.3 
##  ABG_5.563  ABG_5.608  ABG_5.656  ABG_5.674  ABG_5.879  ABG_5.880  ABG_5.881 
##       14.4        5.2        5.7       10.0       10.4        6.7        8.5 
##  ABG_5.882  ABG_5.883  ABG_5.884  ABG_5.887  ABG_5.888  ABG_5.899  ABG_5.900 
##        7.9       11.5        9.0       11.0       10.0       10.3        5.1
```

`obsRecGrowth()` calculates the growth from stems which recruited within
a census interval. 

There are multiple methods for estimating the growth of recruits. The default
`rec_method = "zero"`, assumes that the stem grew from a size of zero in the
initial census. While this method is very simple, it is often unrealistic in
instances where the minimum diameter threshold above which stems are measured
is greater than zero. 


```r
obsRecGrowth(bicuar_clean, "2019", "2021", w = "diam", group = c("plot_id",
    "stem_id"), census = "census_date", rec_method = "zero")
```

```
## ABG_15.929 ABG_15.930  ABG_5.679  ABG_5.680  ABG_5.681  ABG_5.682  ABG_5.683 
##        9.2        8.4        5.3        5.2        5.0        5.2        5.4 
##  ABG_5.684  ABG_5.685  ABG_5.686  ABG_5.687  ABG_5.688 
##        6.2        6.3        5.0        5.0        6.2
```

`rec_method = "thresh"` assumes that the stem grew from the minimum diameter
threshold of the plot. This method works well when growth is being calculated
in terms of diameter, as the value of `w_min_diam` can simply be a numeric
value defining the minimum diameter threshold of the plot. 


```r
obsRecGrowth(bicuar_clean, "2019", "2021", w = "diam", group = c("plot_id",
    "stem_id"), census = "census_date", rec_method = "thresh",
    w_min_diam = 5)
```

```
## ABG_15.929 ABG_15.930  ABG_5.679  ABG_5.680  ABG_5.681  ABG_5.682  ABG_5.683 
##        4.2        3.4        0.3        0.2        0.0        0.2        0.4 
##  ABG_5.684  ABG_5.685  ABG_5.686  ABG_5.687  ABG_5.688 
##        1.2        1.3        0.0        0.0        1.2
```

When `w` is something other than diameter, biomass for instance, this method
requires extra work. In this case, it is necessary to estimate what the biomass
of each stem would be if it shrunk down to the diameter of the minimum diameter
threshold of the plot. In the example below, biomass is estimated using a
simplistic allometric equation.


```r
# Define minimum diameter
mdt <- 5

# Estimate AGB at minimum diameter
bicuar_clean$agb_min <- (0.0644 * (pi * (mdt/2)^2)^1.3341)/1000

# Estimate growth from recruits
obsRecGrowth(bicuar_clean, "2019", "2021", w = "diam", group = c("plot_id",
    "stem_id"), census = "census_date", rec_method = "thresh",
    w_min_diam = "agb_min")
```

```
## ABG_15.929 ABG_15.930  ABG_5.679  ABG_5.680  ABG_5.681  ABG_5.682  ABG_5.683 
##      9.197      8.397      5.297      5.197      4.997      5.197      5.397 
##  ABG_5.684  ABG_5.685  ABG_5.686  ABG_5.687  ABG_5.688 
##      6.197      6.297      4.997      4.997      6.197
```

Finally, `rec_method = "extrap"` uses the average growth rate of the minimum
diameter size class to extrapolate the size of the recruit back to its
estimated size in the initial census. This method relies on having enough data
of sufficient quality to calculate an accurate 'average' growth rate for small
stems. `growth_percentile` refers to the percentile of growth rates used to
calculate the average growth rate of small stems. The default of 0.86 is
recommended by [@Talbot2014].


```r
obsRecGrowth(bicuar_clean, "2019", "2021", w = "diam", group = c("plot_id",
    "stem_id"), census = "census_date", rec_method = "extrap",
    diam = "diam", min_size_class = c(5, 10), min_diam_thresh = 5,
    growth_percentile = 0.86)
```

```
## ABG_15.929 ABG_15.930  ABG_5.684  ABG_5.685  ABG_5.688 
##        9.2        8.4        6.2        6.3        6.2
```

# Average growth rates

Rather than estimating individual stem growth rates separately for each census
interval, as in `obsSurGrowth()`, you may want to calculate the average growth
rate for each individual across all available census intervals. `growthMod()`
fits a linear model to each stem measurement timeline and extracts the model
slope coefficient, which is the average stem growth rate, in units of `w` per
`census`. Like `obsSurGrowth()`, with option `full = TRUE`, `growthMod()`
returns a dataframe with the columns in `group`, rather than a named vector.


```r
growthMod(bicuar_clean, w = "diam", group = c("plot_id", "stem_id"),
    census = "census_date")
```

```
##    ABG_5.1   ABG_5.10  ABG_5.100 ABG_15.101  ABG_5.101  ABG_5.102 ABG_15.103 
##  3.336e-01  8.878e-01  2.147e-01  2.753e-01  1.008e+00  5.544e-01 -1.000e-01 
##  ABG_5.103 ABG_15.104  ABG_5.104 ABG_15.105  ABG_5.105 ABG_15.106  ABG_5.106 
## -4.000e-01 -3.500e-01  1.720e-01  6.213e-01  5.787e-01  1.000e-01  6.458e-01 
## ABG_15.107  ABG_5.107 ABG_15.108  ABG_5.108 ABG_15.109  ABG_5.109   ABG_5.11 
##  8.742e-01  9.363e-01  9.426e-01  1.091e+00  5.350e+00 -1.450e+00  4.844e-01 
## ABG_15.110  ABG_5.110 ABG_15.111  ABG_5.111 ABG_15.112  ABG_5.112 ABG_15.113 
##  1.049e-01  2.341e-01  2.000e-01  3.866e-01  1.090e+00  1.087e+00  1.045e+00 
##  ABG_5.113 ABG_15.114 ABG_15.115  ABG_5.115 ABG_15.116  ABG_5.116 ABG_15.117 
##  2.075e-01  8.594e-01 -5.000e-02  2.517e-01  5.024e-15  6.169e-01  2.512e-15 
##  ABG_5.117 ABG_15.118  ABG_5.118 ABG_15.119  ABG_5.119   ABG_5.12 ABG_15.120 
##  2.996e-01 -1.000e-01  6.373e-01  2.655e-01  1.498e+00  1.596e-01  2.100e+00 
##  ABG_5.120  ABG_5.121 ABG_15.122  ABG_5.122 ABG_15.123  ABG_5.123 ABG_15.124 
## -1.000e-01  2.075e-01  4.000e-01  7.342e-01  1.008e+00  6.245e-01  1.534e-01 
##  ABG_5.124 ABG_15.125  ABG_5.125  ABG_5.126 ABG_15.127  ABG_5.127 ABG_15.128 
##  1.677e-01  7.113e-01  2.482e-01  3.191e-01  8.981e-01  6.500e-01  0.000e+00 
##  ABG_5.128 ABG_15.129  ABG_5.129 ABG_15.130  ABG_5.130 ABG_15.131  ABG_5.131 
##  1.670e-01  1.021e+00  7.016e-01  2.500e-01  1.272e+00  7.763e-01  3.040e-01 
## ABG_15.132  ABG_5.132  ABG_5.133 ABG_15.134  ABG_5.134 ABG_15.135  ABG_5.135 
##  5.000e-02  9.552e-01  6.277e-01  7.402e-01  4.770e-01 -1.000e-01  8.543e-01 
## ABG_15.136  ABG_5.136  ABG_5.137 ABG_15.138  ABG_5.138 ABG_15.139  ABG_5.139 
##  9.637e-01  1.111e+00  9.087e-01  7.544e-01  1.200e+00  2.500e-01  1.180e+00 
##   ABG_5.14 ABG_15.140  ABG_5.140 ABG_15.141  ABG_5.141 ABG_15.142  ABG_5.142 
##  3.536e-01 -5.000e-02  1.314e+00  2.512e-15  1.114e+00  3.994e-01  1.338e+00 
## ABG_15.143  ABG_5.143 ABG_15.144  ABG_5.144 ABG_15.145  ABG_5.145 ABG_15.146 
##  8.913e-01  3.155e-01  0.000e+00  7.086e-01 -1.000e-01  1.005e+00  2.000e-01 
##  ABG_5.146 ABG_15.147  ABG_5.147 ABG_15.148  ABG_5.148 ABG_15.149  ABG_5.149 
##  1.004e+00  5.024e-15  5.000e-01  2.981e-01  1.130e+00  5.000e-02  1.736e+00 
##   ABG_5.15 ABG_15.150 
##  4.625e-01  1.500e-01 
##  [ reached getOption("max.print") -- omitted 703 entries ]
```

With the `growthAll()` function, you can also calculate growth increment for
every available consecutive or pairwise census interval for each stem. Where a
stem measurement in a particular census is missing, the function simply
calculates the growth rate using the next non-missing value:


```r
growthAll(bicuar_clean, w = "diam", group = c("plot_id", "stem_id"),
    census = "census_date", type = "pairwise")
```

```
##                plot_id stem_id   t0   tT int   w0     wT      g
## ABG_5.1.1        ABG_5       1 2019 2021   2 14.6 15.000 0.4000
## ABG_5.1.2        ABG_5       1 2019 2023   4 14.6 15.934 1.3345
## ABG_5.1.3        ABG_5       1 2021 2023   2 15.0 15.934 0.9345
## ABG_5.10.4       ABG_5      10 2019 2021   2 22.2 23.200 1.0000
## ABG_5.10.41      ABG_5      10 2019 2023   4 22.2 25.751 3.5510
## ABG_5.10.42      ABG_5      10 2021 2023   2 23.2 25.751 2.5510
## ABG_5.100.7      ABG_5     100 2019 2021   2  7.0  7.000 0.0000
## ABG_5.100.71     ABG_5     100 2019 2023   4  7.0  7.859 0.8587
## ABG_5.100.72     ABG_5     100 2021 2023   2  7.0  7.859 0.8587
## ABG_15.101.10   ABG_15     101 2019 2021   2 28.9 29.100 0.2000
## ABG_15.101.101  ABG_15     101 2019 2023   4 28.9 30.001 1.1011
## ABG_15.101.102  ABG_15     101 2021 2023   2 29.1 30.001 0.9011
##  [ reached 'max' / getOption("max.print") -- omitted 2051 rows ]
```


```r
growthAll(bicuar_clean, w = "diam", group = c("plot_id", "stem_id"),
    census = "census_date", type = "consecutive")
```

```
##                plot_id stem_id   t0   tT int   w0     wT       g
## ABG_5.1.1        ABG_5       1 2019 2021   2 14.6 15.000  0.4000
## ABG_5.1.2        ABG_5       1 2021 2023   2 15.0 15.934  0.9345
## ABG_5.10.4       ABG_5      10 2019 2021   2 22.2 23.200  1.0000
## ABG_5.10.41      ABG_5      10 2021 2023   2 23.2 25.751  2.5510
## ABG_5.100.7      ABG_5     100 2019 2021   2  7.0  7.000  0.0000
## ABG_5.100.71     ABG_5     100 2021 2023   2  7.0  7.859  0.8587
## ABG_15.101.10   ABG_15     101 2019 2021   2 28.9 29.100  0.2000
## ABG_15.101.101  ABG_15     101 2021 2023   2 29.1 30.001  0.9011
## ABG_5.101.12     ABG_5     101 2019 2021   2  5.2  5.100 -0.1000
## ABG_5.101.121    ABG_5     101 2021 2023   2  5.1  9.230  4.1304
## ABG_5.102.16     ABG_5     102 2019 2021   2 57.3 57.200 -0.1000
## ABG_5.102.161    ABG_5     102 2021 2023   2 57.2 59.518  2.3176
##  [ reached 'max' / getOption("max.print") -- omitted 1421 rows ]
```

# Areal growth rates

Measuring stem growth per unit area is best accomplished using the basal area
of all stems within the plot. First, calculate the basal area of each stem
measurement in units of m^2^, then aggregate to the total basal area within
each plot, and finally use the `growthAll()` function to calculate the basal
area growth for each census interval.


```r
bicuar_plot_ba <- bicuar_clean %>%
    mutate(ba = (pi * (diam/2)^2)/10000) %>%
    group_by(plot_id, census_date) %>%
    summarise(ba_sum = sum(ba, na.rm = TRUE))
```

```
## `summarise()` has grouped output by 'plot_id'. You can override using the
## `.groups` argument.
```

```r
growthAll(bicuar_plot_ba, w = "ba_sum", group = "plot_id", census = "census_date",
    type = "consecutive")
```

```
##          plot_id   t0   tT int     w0     wT        g
## ABG_15.1  ABG_15 2019 2021   2 10.637 10.702  0.06507
## ABG_15.2  ABG_15 2021 2023   2 10.702  6.152 -4.54948
## ABG_5.4    ABG_5 2019 2021   2  7.897  7.125 -0.77190
## ABG_5.41   ABG_5 2021 2023   2  7.125  8.676  1.55031
```

# Estimating unseen growth 

To accurately estimate productivity within a given census interval, i.e. the
total amount of biomass produced, as well as the observed growth of survivors
and the observed growth of recruits it is important to factor in the growth
from stems which died some time between the two censuses. [@Talbot2014]
provides a method for estimating this growth, assuming that the stems die
halfway through the census interval, and that they grow at the median growth
rate of their diameter size class.


```r
unobsMorGrowth(bicuar_clean, "2019", "2021", w = "diam", group = c("plot_id",
    "stem_id"), census = "census_date", diam = "diam")
```

```
##  ABG_5.328  ABG_5.879  ABG_5.360  ABG_5.422  ABG_5.674  ABG_5.883  ABG_5.406 
##      0.000      0.000      0.000      0.000      0.000      0.000      0.000 
##  ABG_5.887  ABG_5.888  ABG_5.899  ABG_5.563 ABG_15.121 ABG_15.126  ABG_5.339 
##      0.000      0.000      0.000      0.000      0.025      0.000      0.000 
##   ABG_5.13  ABG_5.345  ABG_5.226  ABG_5.228  ABG_5.239   ABG_5.24  ABG_5.256 
##      0.000      0.000      0.000      0.000      0.000      0.000      0.000 
##  ABG_5.261  ABG_5.322  ABG_5.326  ABG_5.393  ABG_5.335  ABG_5.881  ABG_5.343 
##      0.000      0.000      0.000      0.000      0.000      0.000      0.000 
##  ABG_5.464  ABG_5.347   ABG_5.49  ABG_5.361  ABG_5.362  ABG_5.363  ABG_5.366 
##      0.000      0.000      0.000      0.000      0.000      0.000      0.000 
##  ABG_5.368  ABG_5.392  ABG_5.549  ABG_5.880  ABG_5.608  ABG_5.453  ABG_5.882 
##      0.000      0.000      0.000      0.000      0.000      0.000      0.000 
##  ABG_5.475  ABG_5.884  ABG_5.540  ABG_5.900  ABG_5.656  ABG_5.203  ABG_5.114 
##      0.000      0.000      0.000      0.000      0.000      0.200      0.200
```

Additionally, there may be some stems which both recruit and then subsequently
die within a census interval, meaning they are never measured. Again,
[@Talbot2014] provides a method to both estimate the number of unseen recruits
which died (`unobsRec()`), and estimate their growth (`unobsRecGrowth()`). This
function only provides a single value, so we must subset the data to a single
plot and single census interval.


```r
bicuar_p1 <- bicuar_clean %>%
    filter(plot_id == "ABG_5")

unobsRec(bicuar_p1, "2019", "2021", group = "stem_id", census = "census_date")
```

```
## [1] 0.3219
```


```r
unobsRecGrowth(bicuar_p1, "2019", "2021", w = "diam", group = "stem_id",
    census = "census_date", diam = "diam", min_size_class = c(5,
        10), w_min_diam = 5)
```

```
## [1] 0
```

Like `obsRecGrowth()` with the argument `rec_method = "extrap"`,
`unobsRecGrowth()` uses the median growth rate of the smallest diameter size
class to estimate the growth rate of the unseen recruits. The function assumes
that these stems recruit in about 1/3 of the way through the census interval,
and die about 2/3 of the way through.

# Plot-level rates of biomass change, productivity, loss

The functions above present parts of a workflow that can be used to generate
estimates of plot-level productivity, biomass turnover, and biomass change.
There are two functions in `{ProdVital}` which provide these plot-level
statistics. `prodKohyama()` uses methods from [@Kohyama2019], while
`prodTalbot()` uses methods from [@Talbot2014]. 

`prodTalbot()` uses `unobsRecGrowth()` and `unobsMorGrowth()` to estimate the
unseen productivity from stems which died, while `prodKohyama()` simply uses
the observed initial and final biomass values of survivors and the population,
and accounts for unseen productivity using estimates of population level
processes.


```r
prodKohyama(bicuar_p1, "2019", "2021", w = "diam", group = "stem_id",
    census = "census_date")
```

```
## $int
## [1] 2
## 
## $t0
## [1] "2019"
## 
## $tT
## [1] "2021"
## 
## $N0
## [1] 691
## 
## $NT
## [1] 654
## 
## $Ns0
## [1] 644
## 
## $Nr0
## [1] 10
## 
## $Nd0
## [1] 47
## 
## $B0
## [1] 6610
## 
## $BT
## [1] 6179
## 
## $Bs0
## [1] 6138
## 
## $Br0
## [1] 0
## 
## $Bd0
## [1] 472
## 
## $dB
## [1] -431.8
## 
## $dB_ann
## [1] -215.9
## 
## $W_max
## [1] 47.12
## 
## $r_turn
## [1] 0.007704
## 
## $m_turn
## [1] 0.03522
## 
## $p_turn
## [1] 0.003265
## 
## $l_turn
## [1] 0.03704
## 
## $Nw
## [1] 672.3
## 
## $Nw_ann
## [1] 681.6
## 
## $Bwk
## [1] 6392
## 
## $Bw_ann
## [1] 6501
## 
## $P_simple
## [1] 20.11
## 
## $L_simple
## [1] 236
## 
## $Bw_simple
## [1] 6394
## 
## $P_ann
## [1] 20.49
## 
## $L_ann
## [1] 236.4
## 
## $P
## [1] 20.87
## 
## $L
## [1] 236.8
## 
## $Bw
## [1] 6392
## 
## $Psimp
## [1] -108143
## 
## $Psimp_clark
## [1] -108143
```


```r
prodTalbot(bicuar_p1, "2019", "2021", w = "agb", diam = "diam",
    group = "stem_id", census = "census_date", w_min_diam = "agb_min",
    rec_method = "zero")
```

```
## $t0
## [1] "2019"
## 
## $tT
## [1] "2021"
## 
## $int
## [1] 2
## 
## $sur_obs
## [1] -0.7859
## 
## $rec_obs
## [1] 0.6452
## 
## $mor_obs
## [1] 7.392
## 
## $mor_unobs
## [1] -0.1136
## 
## $rec_unobs
## [1] 0
## 
## $AGWP_obs
## [1] -0.1407
## 
## $AGWP_est
## [1] -0.2543
## 
## $AGWP_obs_ann
## [1] -0.07034
## 
## $AGWP_est_ann
## [1] -0.1272
```

# Plot-level vital rates

Estimates of population demographic processes are equally important to
understanding ecosystem function as estimates of productivity. `{ProdVital}`
provides `vitalKohyama()` to estimate rates of recruitment, mortality and other
derived metrics, using methods from [@Kohyama2018]. 


```r
vitalKohyama(bicuar_p1, "2019", "2021", group = "stem_id", census = "census_date",
    plot_area = 1)
```

```
## $int
## [1] 2
## 
## $t0
## [1] "2019"
## 
## $tT
## [1] "2021"
## 
## $N0
## [1] 691
## 
## $NT
## [1] 654
## 
## $NST
## [1] 644
## 
## $nT
## [1] 10
## 
## $g
## [1] 0.007704
## 
## $m
## [1] 0.03522
## 
## $r
## [1] 0.007704
## 
## $l
## [1] 0.9729
## 
## $ma
## [1] 0.03461
## 
## $ra
## [1] 0.008299
## 
## $raf
## [1] 0.007675
## 
## $ras
## [1] 0.007734
## 
## $raz
## [1] 0.00721
## 
## $M
## [1] 24.34
## 
## $R
## [1] 5.178
## 
## $Ma
## [1] 23.91
## 
## $Ra
## [1] 5.088
## 
## $Ras
## [1] 5.27
```

# Scaling up to multiple plots and multiple census intervals

All the functions in `{ProdVital}` act on a single plot across a single census
interval, but it is trivial to run the functions across groups of plots and
census intervals using `lapply()` or `for` loops.

First, split the dataset into plots:


```r
bicuar_split <- split(bicuar_clean, bicuar_clean$plot_id)
```

Then, for each plot, identify all the possible pairwise combinations of
censuses:


```r
prod_all <- do.call(rbind, lapply(bicuar_split, function(x) {
    # Find all census dates
    census_date_all <- sort(unique(x$census_date))

    # Create all pairwise combinations of censuses
    comb_list_pair <- combn(census_date_all, 2, simplify = FALSE)

    pK_census <- lapply(comb_list_pair, function(y) {
        # Estimate productivity according to Kohyama
        prodKohyama(x, t0 = y[1], tT = y[2], w = "agb", group = "stem_id",
            census = "census_date")
    })

    pT_census <- lapply(comb_list_pair, function(y) {
        # Estimate productivity according to Kohyama
        prodTalbot(x, t0 = y[1], tT = y[2], w = "agb", diam = "diam",
            group = "stem_id", census = "census_date", w_min_diam = "agb_min",
            rec_method = "zero")
    })

    # Create dataframes of metrics
    pK <- data.frame()[seq_along(pK_census), ]
    for (y in names(pK_census[[1]])) {
        pK[[y]] <- unlist(lapply(pK_census, "[[", y))
    }

    pT <- data.frame()[seq_along(pT_census), ]
    for (y in names(pT_census[[1]])) {
        pT[[y]] <- unlist(lapply(pT_census, "[[", y))
    }

    # Join the dataframes
    out <- left_join(pK, pT, by = c("t0", "tT", "int"))

    # Add plot ID
    out$plot_id <- unique(x$plot_id)

    out
}))

prod_all
```

```
##          int   t0   tT  N0  NT Ns0 Nr0 Nd0    B0    BT   Bs0 Br0     Bd0
## ABG_15.1   2 2019 2021 151 151 149   2   2 69.42 71.10 69.10   0  0.3198
## ABG_15.2   4 2019 2023 151  71  69   2  82 69.42 42.41 33.97   0 35.4530
##               dB  dB_ann W_max   r_turn   m_turn  p_turn   l_turn  Nw Nw_ann
## ABG_15.1   1.672  0.8361 2.966 0.006667 0.006667 0.01421 0.002308 151  151.0
## ABG_15.2 -27.016 -6.7540 2.966 0.007143 0.195793 0.05546 0.178680 106  116.3
##            Bwk Bw_ann P_simple L_simple Bw_simple  P_ann L_ann      P      L
## ABG_15.1 70.26  69.84    0.996   0.1599     70.26 0.9972 0.161 0.9983 0.1622
## ABG_15.2 54.81  58.26    2.109   8.8633     55.92 2.7785 9.532 3.0397 9.7937
##             Bw  Psimp Psimp_clark sur_obs rec_obs mor_obs mor_unobs rec_unobs
## ABG_15.1 70.26  195.7       195.7   1.803 0.18891  0.3198   0.01429         0
## ABG_15.2 54.81 -444.8      -444.8   8.339 0.09849 35.4530   4.05122         0
##          AGWP_obs AGWP_est AGWP_obs_ann AGWP_est_ann plot_id
## ABG_15.1    1.992    2.006        0.996        1.003  ABG_15
## ABG_15.2    8.437   12.488        2.109        3.122  ABG_15
##  [ reached 'max' / getOption("max.print") -- omitted 4 rows ]
```

# Comparing Kohyama and Talbot

Some of the metrics returned by `prodKohyama()` and `prodTalbot()` are
analogous, though their values may differ due to the different methods used by
each function.

For example, comparing annual rates of productivity `P_ann` and `AGWP_est_ann`:


```r
ggplot(prod_all, aes(x = P_ann, y = AGWP_obs_ann)) + geom_abline(linetype = 2,
    colour = "red") + geom_point()
```

<div class="figure" style="text-align: center">
<img src="figure/unnamed-chunk-30-1.png" alt="plot of chunk unnamed-chunk-30"  />
<p class="caption">plot of chunk unnamed-chunk-30</p>
</div>

# References

