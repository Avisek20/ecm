function [arval, nmval, wlival, silval] = max_evals_nsga2(population_size, cent)
    
    global data K dim labels;
    
    clust = zeros(population_size, size(data,1));
    nmval = zeros(1,population_size);
    arval = zeros(1,population_size);
    wlival = zeros(1,population_size);
    silval = zeros(1,population_size);
    
    for cont = 1:population_size
        
        v = reshape(cent(cont,:), dim/K, K)';
        u = zeros(size(data,1),K);
        
        for i = 1:size(data,1)
            
            sum1 = 0;
            for j = 1:K
                sum1 = sum1 + (1/norm(data(i,:)-v(j,:)));
            end
            for j = 1:K
                u(i,j) = (1/norm(data(i,:)-v(j,:))) / sum1;
            end
            
            [~,ind]= max(u(i,:));
            
            clust(cont,i)= ind;
        
        end
        
        nmval(cont)= nmi(labels,clust(cont,:));
        arval(cont)= adjrand(labels,clust(cont,:));
        wlival(cont) = WLI(data, v, u);
        sil_tmp = silhouette(data, clust(cont,:));
        silval(cont) = mean(sil_tmp);
    
    end
    
end