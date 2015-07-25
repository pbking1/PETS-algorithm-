- Run drugScore.m. This file will call computePETOneDrug.
- computePETOneDrug is the main engine and must be optimized.

I) INPUTS
- adj.mat: matrix for protein-protein interaction. Possible values: -1, 1, 0
- proteinName.mat: list of protein names
- drugName.mat: list of drug names
- drugVector.mat: matrix of drug-target. Each row stands for a drug, each column stands for a protein
- Rp.mat: Rp scores for all proteins
- ExpressionER+.mat: expression for all proteins

II) Output
After you run drugScore.m, three will be a variable named drugRankingScore in the Matlab workspace. This is how the output looks like.
