# DirOutClassifier
Non-parametric classifier for functional data based on directional outlyingness.

The method first does a dimensionality reduction using the directional outlyingness from Dai and Genton https://arxiv.org/abs/1612.04615. 
You can choose: 
- 3 different univariate outlyingness measures (Mahalanobis, Stahel-Donoho and Adjusted outlyingness), 
- using the outlyingness measure or convert it to a depth measure,
- choose the L-norm used to calculate the variation of directions.

Secondly, after the dimensionality reduction is done, the method classifies the observations using the DD-plot, where you can choose between 3 different classifiers for multivariate data (KNN, SVM and QDA).

