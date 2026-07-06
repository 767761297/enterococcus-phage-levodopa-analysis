[[1]]
[1] "ANOVA(betadisper):"

[[2]]
Analysis of Variance Table

Response: Distances
          Df Sum Sq Mean Sq F value Pr(>F)
Groups     3  10.84  3.6133  0.3891 0.7624
Residuals 16 148.59  9.2866               

[[3]]
[1] ""

[[4]]
[1] "permutest(betadisper):"

[[5]]

Permutation test for homogeneity of multivariate
dispersions
Permutation: free
Number of permutations: 9999

Response: Distances
          Df Sum Sq Mean Sq      F N.Perm Pr(>F)
Groups     3  10.84  3.6133 0.3891   9999 0.7632
Residuals 16 148.59  9.2866                     

[[6]]
[1] ""

[[7]]
[1] "TukeyHSD(betadisper):"

[[8]]
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = distances ~ group, data = df)

$group
                                diff       lwr      upr
MR-GR                     -0.5074991 -6.021647 5.006649
MR+Phage cocktail-GR      -0.6956583 -6.209806 4.818489
MR+AFMT-GR                 1.1943548 -4.319793 6.708503
MR+Phage cocktail-MR      -0.1881592 -5.702307 5.325989
MR+AFMT-MR                 1.7018540 -3.812294 7.216002
MR+AFMT-MR+Phage cocktail  1.8900131 -3.624135 7.404161
                              p adj
MR-GR                     0.9933481
MR+Phage cocktail-GR      0.9833191
MR+AFMT-GR                0.9242768
MR+Phage cocktail-MR      0.9996519
MR+AFMT-MR                0.8136016
MR+AFMT-MR+Phage cocktail 0.7622360


