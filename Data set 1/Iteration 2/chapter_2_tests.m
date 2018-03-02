old_info = 'model_parameters_free_t1.mat'; 
load(old_info)

% load off target data 
[original_sites, off_target_sites, scores] = load_off_target_data; 
% load primary data 
[average_binding_data, average_sgRNA_Tarasava, average_host_Tarasava] = ...
            load_avg_data;

% re-save off target_stuff
save(old_info,'model_parameters_free')
%% initialization of parameters
load(old_info)
[r_crRNA, Delta_crRNA, r_Cas9, Delta_Cas9, Delta_complex, k_I, k_f, Lambda, D, V, k_d, k_c, mu] = initialize_parameters();
ODE_parameters = [r_crRNA, Delta_crRNA, r_Cas9, Delta_Cas9, Delta_complex, k_I, k_f, Lambda, D, V, k_d, k_c, mu]';
integration_time = (600); % only consider the first 10 minutes  
integration_accuracy = 4; % mutiple of 1 evaluation per minute
Cas_concentration = 25; % concentration of Cas molecule injected to start nanoleters
crRNA_concentration = 25; % concentration of crRNA molecule injected to start nanoleters
on_target_copy_number = 1; % number of duplicates of a host site (off target)
off_target_copy_number = ones(length(scores)+1,1); % number of duplicates of a host site (off target)
percent_off_target = 1; % written as percentage 

% train model parameters  
max_iterations = .1*10^4;
Number_of_Nucleotides = 32;
use_old = 1;
off_target_weight = .00; % the error from off targeting will be worth 1/10 of the true data

% re-save off target_stuff
save(old_info,'model_parameters_free')

%% test values fig 2
add_plot = 0;
load(old_info)
% mu
res = 100;
mu_plotting = zeros(res,2);
values = exp(linspace(-15, 10^(1),res));
for ii = 1:res
    model_parameters_free(2*33+23+12) = values(ii);
    [error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava(8,:), average_sgRNA_Tarasava(8,:),...
        average_binding_data(8,:), model_parameters_free,...
        Number_of_Nucleotides, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, 0, percent_off_target, scores);
    mu_plotting(ii,1) = log(values(ii));
    mu_plotting(ii,2) = P__binding_rate;
end
figure(1)
plot(mu_plotting(:,1),mu_plotting(:,2))
title('On-Target Binding Rate vs Cellular Growth Rate \mu')
xlabel('\mu')
ylabel('On-Target Binding Rate')
hold on
%V
load(old_info)
res = 100;
V_plotting = zeros(res,2);
values = exp(linspace(-10^(1),10^(1),res));
for ii = 1:res
    model_parameters_free(2*33+23+9) = values(ii);
    [error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava(8,:), average_sgRNA_Tarasava(8,:),...
        average_binding_data(8,:), model_parameters_free,...
        Number_of_Nucleotides, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, 0, percent_off_target, scores);
    V_plotting(ii,1) = log(values(ii));
    V_plotting(ii,2) = P__binding_rate;
end
figure(2)
plot(V_plotting(:,1),V_plotting(:,2))
title('On-Target Binding Rate vs Cell Volume V')
xlabel('V')
ylabel('On-Target Binding Rate')
hold on

%% test values fig 3
add_plot = 0;
load(old_info)
% rate of expression
res = 100;
r_plotting = zeros(res,2);
values = exp(linspace(-10^(1), 10^(1),res));
for ii = 1:res
    model_parameters_free(2*33+23+0) = values(ii);
    model_parameters_free(2*33+23+0) = values(ii);
    [error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava(8,:), average_sgRNA_Tarasava(8,:),...
        average_binding_data(8,:), model_parameters_free,...
        Number_of_Nucleotides, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, 0, percent_off_target, scores);
    r_plotting(ii,1) = log(values(ii));
    r_plotting(ii,2) = P__binding_rate;
end
figure(1)
plot(r_plotting(:,1),r_plotting(:,2))
title('On-Target Binding Rate vs gRNA and dCas9 Expression Rate')
xlabel('r_{dCas9},r_{crRNA}')
ylabel('On-Target Binding Rate')
hold on
% rate of deterioration
load(old_info)
res = 100;
delta_plotting = zeros(res,2);
values = exp(linspace(-10^(1),10^(1),res));
for ii = 1:res
    model_parameters_free(2*33+23+1) = values(ii);
    model_parameters_free(2*33+23+3) = values(ii);
    model_parameters_free(2*33+23+4) = values(ii);
    [error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava(8,:), average_sgRNA_Tarasava(8,:),...
        average_binding_data(8,:), model_parameters_free,...
        Number_of_Nucleotides, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, 0, percent_off_target, scores);
    delta_plotting(ii,1) = log(values(ii));
    delta_plotting(ii,2) = P__binding_rate;
end
figure(2)
plot(delta_plotting(:,1),delta_plotting(:,2))
title('On-Target Binding Rate vs gRNA and dCas9 Deterioration Rate')
xlabel('\delta_{dCas9},\delta_{crRNA}')
ylabel('On-Target Binding Rate')
hold on

%% test values fig 4
add_plot = 0;
load(old_info)
% k_i
res = 100;
k_i_plotting = zeros(res,2);
values = exp(linspace(-10^(1),10^(1),res));
for ii = 1:res
    model_parameters_free(2*33+23+5) = values(ii);
    [error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava(8,:), average_sgRNA_Tarasava(8,:),...
        average_binding_data(8,:), model_parameters_free,...
        Number_of_Nucleotides, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, 0, percent_off_target, scores);
    k_i_plotting(ii,1) = log(values(ii));
    k_i_plotting(ii,2) = P__binding_rate;
end
figure(1)
plot(k_i_plotting(:,1),k_i_plotting(:,2))
title('On-Target Binding Rate vs Coefficient of Isomerization')
xlabel('k_i')
ylabel('On-Target Binding Rate')
hold on
load(old_info)
% k_d
res = 100;
k_i_plotting = zeros(res,2);
values = exp(linspace(-10^(1),10^(1),res));
for ii = 1:res
    model_parameters_free(2*33+23+6) = values(ii);
    [error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava(8,:), average_sgRNA_Tarasava(8,:),...
        average_binding_data(8,:), model_parameters_free,...
        Number_of_Nucleotides, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, 0, percent_off_target, scores);
    k_i_plotting(ii,1) = log(values(ii));
    k_i_plotting(ii,2) = P__binding_rate;
end
figure(2)
plot(k_i_plotting(:,1),k_i_plotting(:,2))
title('On-Target Binding Rate vs Coefficient of Formation')
xlabel('k_d')
ylabel('On-Target Binding Rate')
hold on

%% test values fig 4
add_plot = 0;
load(old_info)
close all
% concentrations
res = 250;
cencentration_plotting = zeros(res,2);
values = linspace(1,2500,res);
for ii = 1:res
    Cas_concentration = values(ii);
    crRNA_concentration = values(ii);
    [error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava(8,:), average_sgRNA_Tarasava(8,:),...
        average_binding_data(8,:), model_parameters_free,...
        Number_of_Nucleotides, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, 0, percent_off_target, scores);
    cencentration_plotting(ii,1) = log(values(ii));
    cencentration_plotting(ii,2) = P__binding_rate;
end
figure(1)
plot(cencentration_plotting(:,1),cencentration_plotting(:,2))
title('On-Target Binding Rate vs Equimolor Concentration')
xlabel('Equimolor Concentration')
ylabel('On-Target Binding Rate')
hold on

add_plot = 0;
load(old_info)
% copy number
res = 5;
copy_number_plotting = zeros(res,2);
values = linspace(1,5,res);
for ii = 1:res
    on_target_copy_number = values(ii);
    [error, P_on_target, P__binding_rate] = test_model(average_host_Tarasava(8,:), average_sgRNA_Tarasava(8,:),...
        average_binding_data(8,:), model_parameters_free,...
        Number_of_Nucleotides, ...
        off_target_sites, integration_time, integration_accuracy,...
        Cas_concentration, crRNA_concentration, on_target_copy_number, ...
        off_target_copy_number, 0, percent_off_target, scores);
    copy_number_plotting(ii,1) = log(values(ii));
    copy_number_plotting(ii,2) = P__binding_rate;
end
 
figure(2)
plot(cencentration_plotting(:,1),cencentration_plotting(:,2))
title('On-Target Binding Rate vs On-Target Copy Number')
xlabel('On Target Copy Number')
ylabel('On-Target Binding Rate')
hold on