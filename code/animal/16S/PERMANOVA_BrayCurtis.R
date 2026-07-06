Permutation test for adonis under reduced model
Permutation: free
Number of permutations: 9999

vegan::adonis2(formula = as.matrix(dist_bc) ~ group, data = df_meta_for_vegan, permutations = 9999)
         Df SumOfSqs     R2      F Pr(>F)    
Model     3   1.6761 0.4349 4.1046  2e-04 ***
Residual 16   2.1779 0.5651                  
Total    19   3.8540 1.0000                  
---
Signif. codes:  0 ¡®***¡¯ 0.001 ¡®**¡¯ 0.01 ¡®*¡¯ 0.05 ¡®.¡¯ 0.1 ¡® ¡¯ 1
