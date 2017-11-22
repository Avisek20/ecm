function [cent] = get_centers_nsga2(population_size)
    global chromosome feat K ;
    disp(size(chromosome));
    % Reconstruct centers from chromosome
    cent = zeros(population_size, feat*K);
    for i = 1:population_size
        for j = 1:feat*K
            cent(i,j) = chromosome(i,j);
        end
    end
end