# hmmSSF

R package to fit state-switching step selection functions (HMM-SSFs) with covariate-dependent transition probabilities. 

The HMM-SSF is a hidden Markov model, where the (behavioural) states are defined by a movement and habitat selection model (i.e., an SSF). In `hmmSSF`, the model is implemented using importance sampling and direct numerical optimisation via the forward algorithm. The package also contains tools to plot results, and conduct local and global decoding of the state process.

Get started with the main vignette: 

[Fit state-switching step selection functions in hmmSSF](https://github.com/NJKlappstein/hmmSSF/blob/main/vignettes/hmmSSF_introduction.pdf)

### Installation
You can install the package with `devtools`:
```{r}
devtools::install_git("NJKlappstein/hmmSSF")
```


### References

Full model and implementation details found in:  

[Klappstein NJ, L Thomas, & T Michelot. 2023. Flexible hidden Markov models for behaviour-dependent habitat selection. Movement Ecology 11:30 ](https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-023-00392-3).

Also, see the paper where the HMM-SSF was first described: 

[Nicosia A, T Duchesne, L-P Rivest, & D Fortin. 2017. Annals of Applied Statistics 11: 1537-1560](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-11/issue-3/A-multi-state-conditional-logistic-regression-model-for-the-analysis/10.1214/17-AOAS1045.full).



