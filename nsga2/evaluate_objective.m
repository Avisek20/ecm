function f = evaluate_objective(x, M, V)
global K dim

z = reshape(x(1:dim),dim/K,K)';

[f(1), f(2)] = entropy_kmeans_mod(z);

%% Check for error
if length(f) ~= M
    error('The number of decision variables does not match you previous input. Kindly check your objective function');
end
