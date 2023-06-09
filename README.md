# hmmSSF

R package to fit state-switching step selection functions (HMM-SSFs) with covariate-dependent transition probabilities. 

The HMM-SSF is a hidden Markov model, where the (behavioural) states are defined by a movement and habitat selection model (i.e., an SSF). In `hmmSSF`, the model is implemented using importance sampling and direct numerical optimisation via the forward algorithm. The package also contains tools to plot results, and conduct local and global decoding of the state process.

### Vignette

The vignette is the best place to get started; it presents an example analysis from data preparation to interpretation and visualisation of results.

- [Fit state-switching step selection functions in hmmSSF](https://github.com/NJKlappstein/hmmSSF/blob/main/vignettes/hmmSSF_introduction.pdf)

### Installation

You can install the package from Github with `devtools`:
``` R
devtools::install_github("NJKlappstein/hmmSSF")
```

Note that the package is under development, so please keep an eye out for future versions (particularly for bug fixes!)

### References

Klappstein NJ, L Thomas, & T Michelot (2023). [Flexible hidden Markov models for behaviour-dependent habitat selection](https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-023-00392-3). Movement Ecology 11:30.

Nicosia, A., Duchesne, T., Rivest, L. P., & Fortin, D. (2017). [A multi-state conditional logistic regression model for the analysis of animal movement](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-11/issue-3/A-multi-state-conditional-logistic-regression-model-for-the-analysis/10.1214/17-AOAS1045.full). Annals of Applied Statistics 11:1537-1560.



