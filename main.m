%% main 

% load off target data 
[original_sites, off_target_sites, scores, copy_number, pam_sites] = load_off_target_data();

% load primary data 
[average_binding_data, average_sgRNA_Tarasava, average_host_Tarasava] = ...
            load_avg_data;
        
% calc sd for each data point
[binding_frequency, Data_matrix_SgRNA, Data_matrix_Host,...
            Data_vector_DataSet, Data_vector_CasType] = load_seq_data;
        
data_instance = [binding_frequency(1:30),binding_frequency(31:60),binding_frequency(61:90)];
sd = ((max(data_instance') - min(data_instance'))./4).^.5;


%% initialization of parameters
[r_crRNA, Delta_crRNA, r_Cas9, Delta_Cas9, Delta_complex, k_I, k_f, Lambda, D, V, k_d, k_c, mu] = initialize_parameters();
ODE_parameters = [r_crRNA, Delta_crRNA, r_Cas9, Delta_Cas9, Delta_complex, k_I, k_f, Lambda, D, V, k_d, k_c, mu]';
integration_time = (10^4); % now consider steady state conditions
integration_accuracy = 1; % mutiple of 1 evaluation per minute
Cas_concentration = 25; % concentration of Cas molecule injected to start nanoleters
crRNA_concentration = 25; % concentration of crRNA molecule injected to start nanoleters
on_target_copy_number = 1; % number of duplicates of a host site (off target)
off_target_copy_number = copy_number; % number of duplicates of a host site (off target)
percent_off_target = 1; % written as percentage 
%[GG, GC, GA, GT, CG, CC, CA, CT, AG, AC, AA, AT, TG, TC, TA, TT ]; 
% second value is the gRNA first is the host 
match_mismatch_params = -1*[-1.97, -2.70, -1.66, -2.04,...
                             -1.44, -1.97, -0.78, -1.29,...
                             -1.29, -2.04, -1.04, -1.27,...
                             -0.78, -1.66, -0.12, -1.04]'; 
% ['ATG', 'GAG', 'GGA', 'GAA', 'GAT']
pam_parameters = -1*[-6.6, -6.9, -9.4, -6.8, -6.9]/max(-1*[-6.6, -6.9, -9.4, -6.8, -6.9]);
% train model parameters  
max_iterations = 2*10^3;
Number_of_Nucleotides = 32;
stopping_Tol = .1;
off_target_weight = .00; % the error from off targeting will be worth 1/10 of the true data
old_info = 'model_parameters_free.mat'; 

%% train 
% call function
use_old = 1;
[model_parameters, just_trained] = train_model(average_host_Tarasava, average_sgRNA_Tarasava,...
                            average_binding_data, max_iterations,...
                            use_old, Number_of_Nucleotides,...
                            old_info, ODE_parameters, off_target_sites,...
                            integration_time,integration_accuracy,...
                            Cas_concentration, crRNA_concentration, ...
                            on_target_copy_number, off_target_copy_number, pam_sites, ...
                            percent_off_target, off_target_weight, stopping_Tol);
%save parameters

model_parameters_free = model_parameters; 
save(old_info, 'model_parameters_free')
save('model_parameters_free1.mat', 'just_trained')
disp(model_parameters_free)


%% test values
add_plot = 0;
load(old_info)
percent_off_target = 1;
[error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava, average_sgRNA_Tarasava,...
        average_binding_data, just_trained, Number_of_Nucleotides, ...
                            off_target_sites, integration_time, integration_accuracy,...
                            Cas_concentration, crRNA_concentration, on_target_copy_number, ...
                            off_target_copy_number, off_target_weight, ODE_parameters, ...
                            pam_sites, add_plot);
    
close all
% accuracy plot on real data
hold on
figure(1)   
b = bar(linspace(1,length(P__binding_rate),length(P__binding_rate)),[P__binding_rate average_binding_data]);
b(2).FaceColor = 'r';
errorbar(1:30, average_binding_data, sd, '-')
title('observed and predicted binding')
xlabel('distance from current position backward')
ylabel('function evaluation')
title('True versus Predicted Bind ing within dCas9 System')
xlabel('Data Point')
ylabel('Predicted Binding')

%% reduce number of parameters via fourier series

current_parameters = just_trained(6:end);
[a,b,yfit] = Fseries(linspace(-1,1,67),current_parameters,12);
yfit(1) = just_trained(1); yfit(end) = just_trained(end);
new_params = [just_trained(1:5), yfit];

add_plot = 0;
load(old_info)
percent_off_target = 1;
[error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava, average_sgRNA_Tarasava,...
        average_binding_data, new_params, Number_of_Nucleotides, ...
                            off_target_sites, integration_time, integration_accuracy,...
                            Cas_concentration, crRNA_concentration, on_target_copy_number, ...
                            off_target_copy_number, off_target_weight, ODE_parameters, ...
                            pam_sites, add_plot);
disp(error)
close all
% accuracy plot on real data
hold on
figure(1)   
b = bar(linspace(1,length(P__binding_rate),length(P__binding_rate)),[P__binding_rate average_binding_data]);
b(2).FaceColor = 'r';
errorbar(1:30, average_binding_data, sd, '-')
title('observed and predicted binding')
xlabel('distance from current position backward')
ylabel('function evaluation')
title('True versus Predicted Bind ing within dCas9 System')
xlabel('Data Point')
ylabel('Predicted Binding')

%% look at parameter sensetivity
h = .1;
gradients = parameter_sensetivity(average_host_Tarasava, average_sgRNA_Tarasava,...
                            average_binding_data, max_iterations,...
                            use_old, Number_of_Nucleotides,...
                            old_info, ODE_parameters, off_target_sites,...
                            integration_time, integration_accuracy,...
                            Cas_concentration, crRNA_concentration, ...
                            on_target_copy_number, off_target_copy_number, pam_sites, ...
                            percent_off_target, off_target_weight, stopping_Tol,h);

%% now do weighted least squares                 
[a0,a,b,yfit] = wls_fourier_approx(linspace(-1,1,67), just_trained(6:end), 40, abs(gradients(6:end,1)));
plot(linspace(-1,1,67), yfit); hold on; plot(linspace(-1,1,67), just_trained(6:end))
new_params = [just_trained(1:5), yfit];
add_plot = 0;
load(old_info)
percent_off_target = 1;
[error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava, average_sgRNA_Tarasava,...
        average_binding_data, new_params, Number_of_Nucleotides, ...
                            off_target_sites, integration_time, integration_accuracy,...
                            Cas_concentration, crRNA_concentration, on_target_copy_number, ...
                            off_target_copy_number, off_target_weight, ODE_parameters, ...
                            pam_sites, add_plot);
disp(error)

                 