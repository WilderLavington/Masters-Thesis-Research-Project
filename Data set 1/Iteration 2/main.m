%% main 

% load off target data 
[original_sites, off_target_sites, scores] = load_off_target_data; 
% load primary data 
[average_binding_data, average_sgRNA_Tarasava, average_host_Tarasava] = ...
            load_avg_data;
        
%% initialization of parameters
[r_crRNA, Delta_crRNA, r_Cas9, Delta_Cas9, Delta_complex, k_I, k_f, Lambda, D, V, k_d, k_c, mu] = initialize_parameters();
ODE_parameters = [r_crRNA, Delta_crRNA, r_Cas9, Delta_Cas9, Delta_complex, k_I, k_f, Lambda, D, V, k_d, k_c, mu]';
integration_time = (600); % only consider the first 10 minutes  
integration_accuracy = 6; % mutiple of 1 evaluation per minute
Cas_concentration = 50; % concentration of Cas molecule injected to start nanoleters
crRNA_concentration = 50; % concentration of crRNA molecule injected to start nanoleters
on_target_copy_number = 1; % number of duplicates of a host site (off target)
off_target_copy_number = ones(length(scores)+1,1); % number of duplicates of a host site (off target)
percent_off_target = .5; % written as percentage 

% train model parameters  
max_iterations = 1*10^4;
Number_of_Nucleotides = 32;
use_old = 1;
off_target_weight = .00; % the error from off targeting will be worth 1/10 of the true data
old_info = 'model_parameters_free_t1.mat'; 

%% train 

% call function
model_parameters_free = train_model(average_host_Tarasava, average_sgRNA_Tarasava,...
                            average_binding_data, max_iterations,...
                            use_old, Number_of_Nucleotides,...
                            old_info, ODE_parameters, off_target_sites,...
                            integration_time,integration_accuracy,...
                            Cas_concentration, crRNA_concentration, ...
                            on_target_copy_number, off_target_copy_number, ...
                            percent_off_target, off_target_weight);
%save parameters
save(old_info, 'model_parameters_free')
disp(model_parameters_free)

%% test values
add_plot = 1;
load(old_info)
percent_off_target = 1;
[error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava, average_sgRNA_Tarasava,...
        average_binding_data, model_parameters_free,...
        Number_of_Nucleotides, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, 1, percent_off_target, scores);
close all
% accuracy plot on real data
figure(1)   
b = bar(linspace(1,length(P__binding_rate),length(P__binding_rate)),[P__binding_rate average_binding_data]);
b(2).FaceColor = 'r';
title('observed and predicted binding')
xlabel('distance from current position backward')
ylabel('function evaluation')
title('True versus Predicted Bind ing within dCas9 System')
xlabel('Data Point')
ylabel('Predicted Binding')



