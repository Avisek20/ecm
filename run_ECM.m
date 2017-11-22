% run_ECM runs ECM-NSGA-II and ECM-MOEA/D on the data sets in the
% 'datasets' folder. The output is written to the 'output' folder. The
% folder 'evaluations' contains the source code for the evaluation indices.
%

% Set constants
folder_name = 'datasets/';
input_paths = 'list_datasets.txt';
output_foler = 'output';
num_iterations = 5000;
tolerance = 1e-16;
population_size = 50;

% Set path to all subfolders
folder = fileparts(which('run_ECM.m')); 
addpath(genpath(folder));

% Read in dataset file names
dataset_file_names = {};
file_names = {};
fr = fopen(input_paths);
tline = fgetl(fr);
iter1 = 1;
while ischar(tline)
    dataset_file_names{iter1} = strcat(strcat(folder_name,tline),'.txt'); %#ok<SAGROW>
    file_names{iter1} = tline; %#ok<SAGROW>
    iter1 = iter1 + 1;
    tline = fgetl(fr);
end
fclose(fr);

global data labels K;
global chromosome dim feat sigma;

fwrite = fopen(strcat(output_foler,'/all_output.txt'),'w');

% Read in datasets one-by-one
for iter1 = 1:size(dataset_file_names,2)
    data = dlmread(dataset_file_names{iter1});
    disp(dataset_file_names{iter1});
    fprintf(fwrite,'\n\nDataset : %s\n',dataset_file_names{iter1});
    
    % Separate out class labels
    labels = data(:,size(data,2)) + 1;
    data = data(:,1:size(data,2)-1);
        
    % max-min normalization
    data = bsxfun(@rdivide, bsxfun(@minus, data, min(data,[],1)), (max(data,[],1)-min(data,[],1)));
    data(isnan(data)) = 0;
    disp(size(data));
    
    % set range to (-1,1)
    data = data*2 -1;
    
    % Compute sigma, to be variance of the data
    sigma = std(sum(bsxfun(@minus,data,mean(data)).^2,2));
    %#if sigma < 0.7
    %    sigma = 0.7;
    %end
    disp('sigma');
    disp(sigma);
    
    
    % Count number of clusters from labels
    K = size(unique(labels),1);
    disp('K'); disp(K);
    fprintf(fwrite,'Size : (%d,%d)\n',size(data,1),size(data,2));
    fprintf(fwrite,'K : %d\n', K);
    
    
    % Run ECM-NSGA-II
    disp('NSGA-II');
    feat = size(data, 2);
    dim = feat*K;
    starttime = clock;
    % Run nsga2
    nsga_2(population_size, num_iterations);
    % Write time taken
    fprintf('Total time spent in execution : %u\n', etime(clock, starttime))
    % Get centers and write to file
    cent_nsga2 = get_centers_nsga2(population_size);
    fprintf(fwrite,'\nECM-NSGA-II\nCenters :\n');
    for iter2=1:size(cent_nsga2,1)
        for iter3=1:size(cent_nsga2,2)
            fprintf(fwrite,'%f ',cent_nsga2(iter2,iter3));
        end
        fprintf(fwrite,'\n');
    end
    % Get pareto front and write to file
    fprintf(fwrite,'Pareto Front :\n');
    pf_nsga2 = get_pf_nsga2(population_size);
    for iter2=1:size(pf_nsga2,1)
        for iter3=1:size(pf_nsga2,2)
            fprintf(fwrite,'%f ',pf_nsga2(iter2,iter3));
        end
        fprintf(fwrite,'\n');
    end
    
    % Evaluate ECM-NSGA-II
    [nsga2_ari, nsga2_nmi, nsga2_wli, nsga2_sil] = max_evals_nsga2(population_size, cent_nsga2);
    % Write ARI to file
    fprintf(fwrite,'ARI : \n');
    for iter2 = 1:population_size
        fprintf(fwrite,'%f ',nsga2_ari(iter2));
    end
    [ari_v,ari_i] = max(nsga2_ari);
    fprintf(fwrite, '\nMax ARI : %f at index : %d\n', ari_v, ari_i);
    % Write NMI to file
    fprintf(fwrite,'NMI : \n');
    for iter2 = 1:population_size
        fprintf(fwrite,'%f ',nsga2_nmi(iter2));
    end
    % Write WLI to file
    fprintf(fwrite,'WLI : \n');
    for iter2 = 1:population_size
        fprintf(fwrite,'%f ',nsga2_wli(iter2));
    end
    [wli_v,wli_i] = min(nsga2_wli);
    fprintf(fwrite, '\nMin WLI : %f at index : %d\n',wli_v,wli_i);
    % Write Silhouette to file
    fprintf(fwrite,'Silhouette : \n');
    for iter2 = 1:population_size
        fprintf(fwrite,'%f ',nsga2_sil(iter2));
    end
    [sil_v,sil_i] = max(nsga2_sil);
    fprintf(fwrite, '\nMax Silhouette : %f at index : %d\n',sil_v,sil_i);
    % Evaluate Schotts Spacing Metric
    spacing_metric_nsga2 = schott(pf_nsga2);
    % Write Spacing Metric to file
    fprintf(fwrite,'Schotts Spacing Metric : %f\n', spacing_metric_nsga2);    
    
    
    % Run ECM-MOEAD
    disp('ECM-MOEAD');
    domain = repmat([min(data); max(data)]', K, 1); % Domain of the decision variables
    ecmMode = 1;
    mop = testmop('ecm', dim, domain, K, data, ecmMode); % Set the mop
    [pareto_moead, pf_moead] = moead( mop, 'popsize', population_size-1, 'niche', 10, ...
        'iteration', num_iterations, 'method', 'te');
    % pf_moead contains pareto front points
    % [pareto_moead.parameter] contains centers
    cent_moead = [pareto_moead.parameter]';
    % Write centers to file
    fprintf(fwrite,'\nECM-MOEAD\nCenters :\n');
    for i = 1:size(cent_moead,1)
        for j = 1:size(cent_moead,2)
            fprintf(fwrite,'%f ',cent_moead(i,j));
        end
        fprintf(fwrite,'\n');
    end
    % Write Pareto Front to file
    fprintf(fwrite,'Pareto Front :\n');
    pf_moead = pf_moead';
    for i=1:size(pf_moead,1)
        for j=1:size(pf_moead,2)
            fprintf(fwrite,'%f ',pf_moead(i,j));
        end
        fprintf(fwrite,'\n');
    end
    % Evaluate ECM-MOEAD
    [moead_ari, moead_nmi, moead_wli, moead_sil] = max_evals_nsga2(population_size, cent_moead);
    % Write ARI to file
    fprintf(fwrite,'ARI : \n');
    for iter2 = 1:population_size
        fprintf(fwrite,'%f ',moead_ari(iter2));
    end
    [ari_v,ari_i] = max(moead_ari);
    fprintf(fwrite, '\nMax ARI : %f at index : %d\n', ari_v, ari_i);
    % Write NMI to file
    fprintf(fwrite,'NMI : \n');
    for iter2 = 1:population_size
        fprintf(fwrite,'%f ',moead_nmi(iter2));
    end
    [nmi_v,nmi_i] = max(moead_nmi);
    fprintf(fwrite, '\nMax NMI : %f at index : %d\n', nmi_v, nmi_i);
    % Write WLI to file
    fprintf(fwrite,'WLI : \n');
    for iter2 = 1:population_size
        fprintf(fwrite,'%f ',moead_wli(iter2));
    end
    [wli_v,wli_i] = min(moead_wli);
    fprintf(fwrite, '\nMin WLI : %f at index : %d\n',wli_v,wli_i);
    % Write Silhouette to file
    fprintf(fwrite,'Silhouette : \n');
    for iter2 = 1:population_size
        fprintf(fwrite,'%f ',moead_sil(iter2));
    end
    [sil_v,sil_i] = max(moead_sil);
    fprintf(fwrite, '\nMax Silhouette : %f at index : %d\n',sil_v,sil_i);
    % Evaluate Schotts Spacing Metric
    spacing_metric_moead = schott(pf_moead);
    % Write Schotts Spacing Metric to file
    fprintf(fwrite,'Schotts Spacing Metric : %f\n', spacing_metric_moead);
    % Compare with ECM-NSGA-II
    [bei,ie12,ie21] = binary_epsilon_indicator(pf_nsga2, pf_moead);
    % Write EIs to file
    fprintf(fwrite,'BEI : %f\nIE(NSGA2,MOEAD) : %f\nIE(MOEAD,NSGA2) : %f\n', bei, ie12, ie21);
    
end

fclose(fwrite);
