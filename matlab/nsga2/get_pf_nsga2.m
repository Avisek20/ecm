function [pf] = get_pf_nsga2(population_size)
    global chromosome feat K ;
    % Reconstruct objective function values
    pf = zeros(population_size, size(chromosome,2)-(feat*K)-2); % Here 2 means the 2 values for tournament selection
    for i = 1:population_size
        for j = (feat*K+1):(size(chromosome,2)-2)
            pf(i,(j-(feat*K)))= chromosome(i,j);
        end
    end
end
