function f = initialize_variables(N, M, V, min_range, max_range)

% K is the total number of array elements. For ease of computation decision
% variables and objective functions are concatenated to form a single
% array. For crossover and mutation only the decision variables are used
% while for selection, only the objective variable are utilized.

f = zeros(N,V+M+2);
%% Initialize each chromosome
% For each chromosome perform the following (N is the population size)
  for i = 1 : N
      for j = 1 : V
          f(i,j) = min_range(j) + (max_range(j) - min_range(j))*rand(1);
      end
    % For ease of computation and handling data the chromosome also has the
    % vlaue of the objective function concatenated at the end. The elements
    % V + 1 to V + M has the objective function valued. 
    % The function evaluate_objective takes one chromosome at a time,
    % infact only the decision variables are passed to the function along
    % with information about the number of objective functions which are
    % processed and returns the value for the objective functions. These
    % values are now stored at the end of the chromosome itself.
    f(i,V+1:V+M) = evaluate_objective(f(i,1:V), M, V);
  end
