function [number_of_objectives, number_of_decision_variables, min_range_of_decision_variable, max_range_of_decision_variable] = objective_description_function()
global K feat data  % initialised in main source file

number_of_objectives = 2;

number_of_decision_variables = feat*K;

% Assumptions : Center components lie in -1 to 1
%min_range_of_decision_variable=-1*ones(1,dim);
% Setting to max and min from dataset instead
%{
min_range_of_decision_variable = zeros(1,feat*K);
tmp1 = min(data,[],1);
for i = 1:K
    start_idx = ((i-1)*feat)+1;
    min_range_of_decision_variable(start_idx:start_idx+feat-1) = tmp1;
end
    
%max_range_of_decision_variable=ones(1,dim);

max_range_of_decision_variable = zeros(1,feat*K);
tmp2 = max(data,[],1);
for i = 1:K
    start_idx = ((i-1)*feat)+1;
    max_range_of_decision_variable(start_idx:start_idx+feat-1) = tmp2;
end
%}

min_range_of_decision_variable = -1*ones(1,feat*K);%min(min(data))*ones(1,feat*K);
max_range_of_decision_variable = ones(1,feat*K);%max(max(data))*ones(1,feat*K);

end
