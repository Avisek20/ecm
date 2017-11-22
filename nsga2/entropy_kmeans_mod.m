function [CK, CE] = entropy_kmeans_mod(z)
global data K sigma

% Initializing some constants
n = size(data,1);
CE = 0; CK = 0;
u = zeros(n,K);
%beta = 10;  % fixing beta at 10

for i = 1:n

    % Computing u
	for j = 1:K
        u(i,j) = exp( -(1/(2*(sigma^2))) * (sum((data(i,:)-z(j,:)).^2)));
	end
    tmp_sum = sum(u(i,:));
    if tmp_sum == 0
        u(i,:) = 1/K;
    else
        u(i,:) = u(i,:)/sum(u(i,:));
    end
   
    % Summing up the objective function values
    for j = 1:K
        CK = CK + u(i,j)*(sum((data(i,:)-z(j,:)).^2));
        if u(i,j) > 1e-16
            CE = CE + u(i,j)*log(u(i,j))/log(2);
        end
    end
    
end


end
