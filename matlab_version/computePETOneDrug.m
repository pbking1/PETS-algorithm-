function [PETScoreMatrix, scoreArray, sMat] = ...
    computePETOneDrug(adj, drugVector, proteinName, Rp, Expression, boostingFactor, sigma, useProteinDegree)
% adj, drugVector, proteinName and drugName should be from the big network

%% the following routine is to compress the RP score
% newMaxRp = 20;
% baseRp = min(Rp);
% oldMaxRp = max(Rp);
% minMult = newMaxRp / oldMaxRp;
% for i = 1:length(Rp)
%     Rp(i,1) = baseRp + minMult * (Rp(i) - baseRp);
% end

%Rp(find(Rp > newMaxRp)) = newMaxRp;

%%
maxIteration = 150;

proteinWeight = Rp;

% create the abs adj matrix
absAdj = abs(adj);

for i = 1:length(proteinName)
    sizeDownStream = length( find(adj(i, :) ~= 0) );
    %sizeUpStream = length( find(adj(:, i) ~= 0) );
    Rp(i) = max(1, sizeDownStream) ;%+ sizeUpStream;
    %Rp(i) = 1 + length( find(adj(i, :) ~= 0) ) ;
    if sizeDownStream > 0
        adj(i, :) = adj(i, :) / sizeDownStream;
    end
    
end

if useProteinDegree == 1
    proteinWeight = Rp;
end



%%
% initialization
numOfPro = length(proteinName);
target = double(drugVector ~= 0);

% now, previous layer and current layer are represented by binary vector
previousLayer = target;
%aaa = absAdj'*previousLayer
currentLayer = sign(absAdj'*previousLayer);
%currentLayer = [];

% initilalize s score, which is drug-protein interaction
scoreArray = zeros(numOfPro, maxIteration);

scoreArray(:, 1) = drugVector .*  Rp;

%%

%keyboard;

%%
% 399 iteration
for k = 2:maxIteration
    %keyboard;
    scoreArray(:, k) = scoreArray(:, k-1);
    
    % from the formula s(i, k) = (1-d)*c + sum( b*d * adj*s(j, k-1) ) if i
    % belongs to the current layer, we try to set up
    
    % the first term: a vector, entry (1-d)*c for component in
    % currentLayer, otherwise 0
    leftTermVector = ((1-sigma)*scoreArray(:, 1)) .* currentLayer;
    % the right term, nodes not belonging to the currentLayer is not
    % updated, therefore create the identity matrix beforehand
    adjThisTurn = eye(numOfPro);
    % now, for the currentLayer row in adjThisTurn, change the row according
    % to adj matrix
    for i = 1:numOfPro
        b_vec = boostingFactor .^ (previousLayer');
        newRow = b_vec .* (sigma*adj(:, i)');
        aa = sigma*adj(:, i)';
        if currentLayer(i) == 1
            adjThisTurn(i, :) = currentLayer(i) * newRow;
            b = adjThisTurn(i, :);
            %1*242 = 1*1 * 1*242
            length(find(b(:,:)~=0))
            %break point in 88 line and i = 20, k = 3 ERROR
        end
    end
    qqq = length(find(adjThisTurn(:,:)~=0));
    %error
    scoreArray(:, k) = leftTermVector + adjThisTurn * scoreArray(:, k-1);
    %242*1 + 242*242 * 242 * 1
    a = leftTermVector + adjThisTurn * scoreArray(:, k-1);
    %keyboard;
    
    % update layers
    previousLayer = currentLayer;
    currentLayer = sign(absAdj'*previousLayer); %error
    
    %PETScoreMatrix(k) = -dot ( scoreArray(:, k),  Expression ) /  sum( abs(scoreArray(:, k)) );
    %keyboard;
    
end
%%

sMat = scoreArray(:, size(scoreArray, 2));
for i = 1:length(proteinName)
    maxSStream = max( abs(scoreArray(i, :)) ); %row
    if abs(scoreArray(i, size(scoreArray, 2)) / maxSStream ) < 0.1 %scoreArray(i, 149)
        sMat(i) = 0;
    end
%     if abs(sMat(i)) > 1
%         sMat(i) = sign(sMat(i));
%     end
end

% compute drug ranking score after the given iteration
nume = 0;
deno = 0;
for i = 1:length(proteinName)
    if Expression(i) ~= 0
        nume = nume - proteinWeight(i) * sign(sMat(i) * Expression(i));
        deno = deno + proteinWeight(i);
    end
end

if deno ~= 0
    PETScoreMatrix = nume / deno;
else
    PETScoreMatrix = 0;
end
PETScoreMatrix
%keyboard;
end