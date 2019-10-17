function subp=init_weights(popsize, niche, objDim)
% init_weights function initialize a pupulation of subproblems structure
% with the generated decomposition weight and the neighbourhood
% relationship.
subp=[];
limit = popsize;
for i=0:limit
    if objDim==2
        p=struct('weight',[],'neighbour',[],'optimal', Inf, 'optpoint',[], 'curpoint', []);
        weight=zeros(2,1);
        weight(1)=i/limit;
        weight(2)=(limit-i)/limit;
        p.weight=weight;
        subp=[subp p];
    elseif objDim==3
        p=struct('weight',[],'neighbour',[],'optimal', Inf, 'optpoint',[], 'curpoint', []);
        weight=zeros(3,1);
        for j = 0:limit
            if (i + j) <= limit
                weight(1) = i / limit;
                weight(2) = j / limit;
                weight(3) = 1 - weight(1) - weight(2);
                p.weight=weight;
                subp=[subp p];
            end
        end
    elseif objDim==4
        p=struct('weight',[],'neighbour',[],'optimal', Inf, 'optpoint',[], 'curpoint', []);
        weight=zeros(4,1);
        for j = 0:limit
            for k = 0:limit
                if (i + j + k) <= limit
                    weight(1) = i / limit;
                    weight(2) = j / limit;
                    weight(3) = k / limit;
                    weight(4) = 1 - weight(1) - weight(2) - weight(3);
                    p.weight=weight;
                    subp=[subp p];
                end
            end
        end
    else
        error('Weight generation method undefined for more than 3 objectives')
    end
    
end

% weight = lhsdesign(popsize, objDim, 'criterion','maximin', 'iterations', 1000)';
% p=struct('weight',[],'neighbour',[],'optimal', Inf, 'optpoint',[], 'curpoint', []);
% subp = repmat(p, popsize, 1);
% cells = num2cell(weight);
% [subp.weight]=cells{:};

%Set up the neighbourhood.
leng=length(subp);
distanceMatrix=zeros(leng, leng);
for i=1:leng
    for j=i+1:leng
        A=subp(i).weight;B=subp(j).weight;
        distanceMatrix(i,j)=(A-B)'*(A-B);
        distanceMatrix(j,i)=distanceMatrix(i,j);
    end
    [s,sindex]=sort(distanceMatrix(i,:));
    subp(i).neighbour=sindex(1:niche)';
end

end