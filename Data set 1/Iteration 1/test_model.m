function [error, P_on_target, P__binding_rate] = test_model(host_site,sg_RNA,binding_rate, ...
            params, Number_of_Nucleotides, ...
            off_target_sites, integration_time, integration_accuracy,...
            Cas_concentration, crRNA_concentration, on_target_copy_number, ...
            off_target_copy_number, add_plot, percent_off_target, scores)
    
    %%%%%========================================
    %set number of off target sites included
    %%%%%========================================
    number_of_off_target_sites = floor(percent_off_target*length(off_target_sites(:,1)));
    off_target_sites = off_target_sites(1:number_of_off_target_sites,:);
    off_target_copy_number = off_target_copy_number(1:number_of_off_target_sites);
    scores = scores(1:number_of_off_target_sites);
    %%%%%========================================
    %set individual model parameters
    %%%%%========================================
    match_mismatch_params = params(1:16);
    anchor_params = params(17:20);
    bias = params(21:22);
    forward = params(23:33+22);
    backword = params(33+23:2*33+22);
    ODE_parameters = params(2*33+23:end);
    
    %%%%%========================================
    %initialization
    %%%%%========================================
    error = 0;
    batch_size = length(binding_rate);
    P_on_target = zeros(batch_size,1);
    P__binding_rate = zeros(batch_size,1);
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
        predicted_binding_rate = 1-x(end,5);
        P__binding_rate(kk) = predicted_binding_rate;
        %%%%%========================================
        %calculate MSE 
        %%%%%========================================
        error = error + ((binding_rate(kk)-predicted_binding_rate)./(mean(binding_rate)))^2;
        
        %%%%%========================================
        %plot values
        %%%%%========================================
        if add_plot
            close all
            hold on
            % plot figures
            figure(1)
        	plot(t, 1-x(:,5))
            title('Percent Binding of On-Target Promotor Site')
            xlabel('Time')
            ylabel('Percentage of Bound Sites')
            figure(2)
            plot(t, x(:,1))
            title('dCas9 Expression')
            xlabel('Time')
            ylabel('Number of Active Molecules')
            figure(3)
            plot(t, x(:,2))
            title('crRNA Expression')
            xlabel('Time')
            ylabel('Number of Active Molecules')
            figure(4)
            plot(t, x(:,3))
            title('Intermediate Complex Expression')
            xlabel('Time')
            ylabel('Number of Active Molecules')
            figure(5)
            plot(t, x(:,4))
            title('Isomerized Complex Expression')
            xlabel('Time')
            ylabel('Number of Active Molecules')
            figure(6)
            for ii = 6:length(x(1,:))
            	plot(t, 1-x(:,ii))
                hold on
            end   
            title('Percent Binding of Off-Target Sites')
            xlabel('Time')
            ylabel('Number of Active Molecules')
            pause
        end
    end
end