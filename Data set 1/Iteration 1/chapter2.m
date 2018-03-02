%% main 

% load off target data 
[original_sites, off_target_sites, scores] = load_off_target_data; 
% load primary data 
[average_binding_data, average_sgRNA_Tarasava, average_host_Tarasava] = ...
            load_avg_data;
        
%% vary equimolor Cas9 gRNA concentration
% ODE info 
[r_crRNA, Delta_crRNA, r_Cas9, Delta_Cas9, Delta_complex, k_I, k_f, Lambda, D, V, k_d, k_c, mu] = initialize_parameters();
ODE_parameters = [r_crRNA, Delta_crRNA, r_Cas9, Delta_Cas9, Delta_complex, k_I, k_f, Lambda, D, V, k_d, k_c, mu]';
integration_time = (600); % only consider the first 10 minutes  
integration_accuracy = 1; % mutiple of 1 evaluation per minute
Cas_concentration = 25; % concentration of Cas molecule injected to start nanoleters
crRNA_concentration = 25; % concentration of crRNA molecule injected to start nanoleters
on_target_copy_number = 1; % number of duplicates of a host site (off target)
off_target_copy_number = ones(size(scores)); % number of duplicates of a host site (off target)
number_of_off_target_sites = 0.1; % written as percentage 

% load old model info 
Number_of_Nucleotides = 32;
fourier_coeffs = 2*(12) + 1; % has to be some odd number 
add_plot = 1;
old_info = 'model_parameters_fourier_localfp25t1.mat';
load(old_info)
% update some of the values 
model_parameters_fourier(2*fourier_coeffs+23:end) = ODE_parameters;
save('model_parameters_fourier_localfp25t1.mat')

%% perfectly-complimentary
[error, P_on_target, predicted_binding_rate] = test_model(average_host_Tarasava(2,:), average_sgRNA_Tarasava(2,:),...
        average_binding_data(2), model_parameters_fourier,...
        Number_of_Nucleotides, fourier_coeffs, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, add_plot,number_of_off_target_sites);
    
%% non-complimentary
test_model(average_host_Tarasava(1,:), average_sgRNA_Tarasava(1,:),...
        average_binding_data(1), model_parameters_fourier,...
        Number_of_Nucleotides, fourier_coeffs, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, add_plot, number_of_off_target_sites);
    
%% partially-complimentary
test_model(average_host_Tarasava(24,:), average_sgRNA_Tarasava(24,:),...
        average_binding_data(24), model_parameters_fourier,...
        Number_of_Nucleotides, fourier_coeffs, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, add_plot,number_of_off_target_sites);
    