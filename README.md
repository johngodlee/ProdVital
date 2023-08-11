# ProdVital

An R package for estimating productivity and vital rates, using repeat measurements of individuals, primarily tree stems, using census-based data from plots. The package draws heavily on work published in Talbot et al. (2014), Kohyama et al. (2018), and Kohyama et al. (2019).

Specifically, the package provides functions to:

1. Identify individuals which survived, recruited or died within a census interval. 
2. Calculate growth increment of survivors and recruits, or losses from individuals which died.
3. Estimate the productivity contribution of stems which died, or recruited and died, within a census interval, using methods from Talbot et al. (2014).
4. Estimate rates of productivity and loss using methods from Talbot et al. (2014) and Kohyama et al. (2019).
5. Estimate rates of individual recruitment and mortality using methods from Kohyama et al. (2018).

## References

Kohyama, Takashi S., Tetsuo I. Kohyama, and Douglas Sheil. 2018. "Definition and Estimation of Vital Rates from Repeated Censuses: Choices, Comparisons and Bias Corrections Focusing on Trees." Edited by Mark Rees. Methods in Ecology and Evolution 9 (4): 809–21. https://doi.org/10.1111/2041-210x.12929.

Kohyama, Takashi S., Tetsuo I. Kohyama, and Douglas Sheil. 2019. "Estimating Net Biomass Production and Loss from Repeated Measurements of Trees in Forests and Woodlands: Formulae, Biases and Recommendations." Forest Ecology and Management 433 (February): 729–40. https://doi.org/10.1016/j.foreco.2018.11.010.

Talbot, Joey, Simon L. Lewis, Gabriela Lopez-Gonzalez, Roel J. W. Brienen, Abel Monteagudo, Timothy R. Baker, Ted R. Feldpausch, et al. 2014. "Methods to Estimate Aboveground Wood Productivity from Long-Term Forest Inventory Plots." Forest Ecology and Management 320: 30–38. https://doi.org/10.1016/j.foreco.2014.02.021.
