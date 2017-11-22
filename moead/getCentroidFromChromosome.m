%% Obtain cluster centroids encoded in chromosomes
function z = getCentroidFromChromosome(x, datasetDim)

K = length(x) / datasetDim; % Number of centroids
z = zeros(K, datasetDim); % Initialize centroid matrix

for centroidNum = 1:K
    
    % Centroids encoded as a series of K vectors [zx, zy] concatenated
    % horizontally
    xIndex = (datasetDim * centroidNum) - (datasetDim - 1);
    z(centroidNum, :) = x(xIndex:xIndex + (datasetDim - 1));
end