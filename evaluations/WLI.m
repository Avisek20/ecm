function [index_wli] = WLI(X, V, U)
	
	k = size(V,1);
	
	WLn = 0;
	for i=1:k
		WLn = WLn + ( (sum((U(:,i).^2).*sum(bsxfun(@minus, X, V(i,:)).^2,2))) / (sum(U(:,i))) );
	end
	
	dist = zeros(1,(k*(k-1))/2);
	dist_idx = 1;
	for i=1:k
		for j=i+1:k
			if i == j
				continue
			end
			dist(dist_idx) = sum((V(i,:) - V(j,:)).^2);
			dist_idx = dist_idx + 1;
		end
	end
	WLd = 0.5*(min(dist) + median(dist));
    
    if ( WLd - 0 ) < 1e-16
        WLd = 1e-16;
    end
	
	index_wli = WLn / (2 * WLd);
end