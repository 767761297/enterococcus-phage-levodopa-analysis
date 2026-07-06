[[1]]
[1] "ANOVA(betadisper):"

[[2]]
Analysis of Variance Table

Response: Distances
          Df Sum Sq Mean Sq F value Pr(>F)
Groups     3 0.8955 0.29849  0.7003 0.5655
Residuals 16 6.8201 0.42626               

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
Groups     3 0.8955 0.29849 0.7003   9999 0.5643
Residuals 16 6.8201 0.42626                     

[[6]]
[1] ""

[[7]]
[1] "TukeyHSD(betadisper):"

[[8]]
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = distances ~ group, data = df)

$group
                                diff        lwr      upr
MR-GR                      0.2279701 -0.9534010 1.409341
MR+Phage cocktail-GR      -0.1742140 -1.3555850 1.007157
MR+AFMT-GR                 0.3789116 -0.8024594 1.560283
MR+Phage cocktail-MR      -0.4021841 -1.5835551 0.779187
MR+AFMT-MR                 0.1509415 -1.0304295 1.332313
MR+AFMT-MR+Phage cocktail  0.5531256 -0.6282454 1.734497
                              p adj
MR-GR                     0.9446290
MR+Phage cocktail-GR      0.9739060
MR+AFMT-GR                0.7958852
MR+Phage cocktail-MR      0.7658578
MR+AFMT-MR                0.9826978
MR+AFMT-MR+Phage cocktail 0.5525780


