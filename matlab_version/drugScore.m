% compute ranking score for each drug

% clear
% clc;

%%
% load the data set

load adj.mat;
load proteinName.mat;
load drugName.mat;
load drugVector.mat;
load Rp.mat;
load ExpressionER+.mat;
%%


%%
% call the computePETOneDrug function for drug ranking
boostingFactor = 1.5; 
% sigma = 0 ; 
 sigma = 0.8;
% user can decide boosing factor and sigma
 
 drugPro = zeros(length(drugName), length(proteinName));
for i = 1:length(drugName)
    % change the last parameter: 0 means using Rp score for weighting, 1
    % mean no Rp score for weighting.
    [PETScoreMatrix, scoreArray, sMat] = ...
        computePETOneDrug(adj, drugVector(:, i), proteinName, Rp, Expression, boostingFactor, sigma, 0);
    drugRankingScore(i, 1) = PETScoreMatrix(length(PETScoreMatrix));
    drugPro(i, :) = sMat';
end

