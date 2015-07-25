###Project description
- this is a project done during the 2015 summer Purdue University MURI program.
- the aim of the project is to convert the matlab version of the program into a Java one
- the original file is PETS and the file in matlab will get the score after running.
- the Java program that I have done is can get exact the same result.
- I use JAMA as the matrix library 
	
####notice
- make sure java is install
- the file of the matlab version should be convert into .csv. Due to the java code can not run a .mat file
	
####input
- the original mat data file is 
	- adj.mat: matrix for protein-protein interaction. Possible values: -1, 1, 0
	- proteinName.mat: list of protein names
	- drugName.mat: list of drug names
	- drugVector.mat: matrix of drug-target. Each row stands for a drug, each column stands for a protein
	- Rp.mat: Rp scores for all proteins
	- ExpressionER+.mat: expression for all proteins
	
	
####output 
- Matlab version 
	- After you run drugScore.m, three will be a variable named drugRankingScore in the Matlab workspace. This is how the output looks like.

- java version of the code will also output the result too and will be the same as the matlab version.

####how to run
- matlab version
	- Run drugScore.m. This file will call computePETOneDrug.
	- computePETOneDrug is the main engine and need be optimized.

- java version	
	- use the following command line
		- sh run.sh
	