function [error] = Obj_function_wo_ODE(binding_rate, sg_RNA, host_site, ...
                            params, Number_of_Nucleotides, ...
                            off_target_sites, integration_time, integration_accuracy,...
                            Cas_concentration, crRNA_concentration, on_target_copy_number, ...
                            off_target_copy_number, off_target_weight, ODE_parameters)
    %%%%%========================================
    %set number of off target sites included
    %%%%%========================================
    number_of_off_target_sites = length(off_target_sites(:,1));

    %%%%%========================================
    %set individual model parameters
    %%%%%========================================
    match_mismatch_params = params(1:16);
    anchor_params = params(17:20);
    bias = params(21:22);
    forward = params(23:33+22);
    backword = params(33+23:2*33+22);
    
    %%%%%========================================
    %initialization
    %%%%%========================================
    error = 0;
    batch_size = length(binding_rate);
    P_on_target = zeros(batch_size,1);

    %%%%%========================================
    %calculate probabilities and set error
    %%%%%========================================
    for kk = 1:batch_size
        
        %%%%%========================================
        % Calculate binding probability for current site
        %%%%%========================================
        N_genome_i = [host_site(kk,:); sg_RNA(kk,:)]';
        [P_lambda, P_mu] = model_probabilities_free(N_genome_i, match_mismatch_params,...
                        anchor_params, Number_of_Nucleotides, ...
                        forward, backword, bias(1), bias(2));
        P_on_target(kk) = BDC_Analytic_Soln(P_lambda,P_mu);
        
        %%%%%========================================
        % Calculate binding probability of off target sites binding
        %%%%%========================================
        batch_size_off_target = length(off_target_sites(:,1));
        P_off_target = zeros(batch_size_off_target,1); 
        for ii = 1:batch_size_off_target 
            N_genome_i = [off_target_sites(ii,:); sg_RNA(kk,:)]';
            [P_lambda, P_mu] = model_probabilities_free(N_genome_i, match_mismatch_params,...
                        anchor_params, Number_of_Nucleotides, ...
                        forward, backword, bias(1), bias(2));
            P_off_target(ii) = BDC_Analytic_Soln(P_lambda,P_mu);
        end
        %%%%%========================================
        % Integrate system of equations to get reletive binding rates
        %%%%%========================================
        Initial_conditions = [Cas_concentration,crRNA_concentration, 0, 0, on_target_copy_number', off_target_copy_number'];
        tspan = linspace(0,integration_time,integration_time*integration_accuracy);
        N_total = [on_target_copy_number', off_target_copy_number'];
        P = [P_on_target(kk), P_off_target'];
        f = @(t,x) Stiff_Solve(t, x, ODE_parameters, N_total, P);
        [t,x] = ode23s(f,tspan,Initial_conditions');
        predicted_binding_rate = on_target_copy_number-x(end,5);
        
        %%%%%========================================
        %Add in scaled off target info for perfectly complimentary site
        %%%%%========================================
        if number_of_off_target_sites ~= 0 
            total_off_target_binding = (sum(off_target_copy_number'-x(end,6:end)));
            rescaled = off_target_weight*(total_off_target_binding/number_of_off_target_sites)^2;
            error = error + rescaled;
        end
        %%%%%========================================
        %calculate MSE 
        %%%%%========================================
        error = error + ((binding_rate(kk)-predicted_binding_rate)./((mean(binding_rate))))^2;
       
    end
    dig = 10^10;
    msg1 = sprintf('%.10f', round(error*dig)/dig);
    msg2 = sprintf('%.2f', length(backword));
    msgs1 = ': ';
    msgs2 = '\n';
    fprintf(msg2);
    fprintf(msgs1);
    fprintf(msg1);
    fprintf(msgs2);
    disp('');
end