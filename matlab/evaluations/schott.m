function [spacing_metric] = schott(pf)
    dist = zeros(size(pf,1),size(pf,1));
    for i = 1:size(pf,1)
        for j = 1:i-1
            dist(i,j) = dist(j,i);
        end
        dist(i,i) = inf;
        for j = i+1:size(pf,1)
            dist(i,j) = sqrt(sum((pf(i,:)-pf(j,:)).^2));
        end
    end
    [mindist, mindist_idx] = min(dist,[],1);
    mean1 = mean(mindist);
    sum1 = 0;
    for i = 1:size(pf,1)
        sum1 = sum1 + sum((mean1 - dist(i, mindist_idx(i))).^2);
    end
    spacing_metric = sqrt(sum1 / (size(pf,1) - 1));
end