function model_parameters = train_model(host_site,sg_RNA,...
                            binding_rate, max_iteration,...
                            use_old, Number_of_Nucleotides,...
                            old_info, ODE_parameters, off_target_sites,...
                            integration_time,integration_accuracy,...
                            Cas_concentration, crRNA_concentration, ...
                            on_target_copy_number, off_target_copy_number, ...
                            percent_off_target, off_target_weight)
    %%%%%========================================
    %Initialize main parameters
    %%%%%========================================
    anchor_params = .5.*rand(4,1);
    bias = .5.*rand(2,1); % [.5;.5]
    free = .5.*rand(2*33,1); 
    ode = ODE_parameters;
   %match_mismatch_params =  2.*rand(16,1);
   %[GG, GC, GA, GT, CG, CC, CA, CT, AG, AC, AA, AT, TG, TC, TA, TT ];
    match_mismatch_params = -1*[-1.97, -2.70, -1.66, -2.04,...
                             -1.44, -1.97, -0.78, -1.29,...
                             -1.29, -2.04, -1.04, -1.27,...
                             -0.78, -1.66, -0.12, -1.04]'; 
    %%%%%========================================
    % check if we are using old info or not
    %%%%%========================================
    if use_old
        load(old_info, 'model_parameters_free');
        current_params = model_parameters_free;
    else
        current_params = [match_mismatch_params; anchor_params; bias; free; ode]';
    end
    
    %%%%%========================================
    %set number of off target sites included random permations
    %%%%%========================================
    number_of_off_target_sites = floor(percent_off_target*length(off_target_sites(:,1)));
    %sites = randperm(length(off_target_sites(:,1)));
    sites = 1:length(off_target_sites(:,1));
    sites = sites(1:number_of_off_target_sites);
    off_target_sites = off_target_sites(sites,:);
    off_target_copy_number = off_target_copy_number(sites);
    
    %%%%%========================================
    %call local optimization algorithm
    %%%%%========================================
    options = optimset('MaxIter', max_iteration,'MaxFunEvals', max_iteration,...
        'TolX', 10^(-64));
    % remove ode and match-mismatch params
    current_params = current_params(17:end);
    current_params = current_params(1:end-13);
    
    % define objective function
    obj_global =  @(params) Obj_function_wo_ODE_MM(binding_rate, sg_RNA, host_site, ...
                            params, Number_of_Nucleotides, ...
                            off_target_sites, integration_time, integration_accuracy,...
                            Cas_concentration, crRNA_concentration, on_target_copy_number, ...
                            off_target_copy_number, off_target_weight, ODE_parameters, match_mismatch_params);
                        
%     obj_global = @(params) Obj_function_wo_ODE(binding_rate, sg_RNA, host_site, ...
%                             params, Number_of_Nucleotides, ...
%                             off_target_sites, integration_time, integration_accuracy,...
%                             Cas_concentration, crRNA_concentration, on_target_copy_number, ...
%                             off_target_copy_number, off_target_weight, ODE_parameters);
    
    % run fminsearch          
    model_parameters = fminsearch(obj_global,current_params,options);
    model_parameters = [match_mismatch_params', model_parameters, ODE_parameters'];
end