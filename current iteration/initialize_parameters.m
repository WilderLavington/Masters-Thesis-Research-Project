function [r_crRNA, Delta_crRNA, r_Cas9, Delta_Cas9, Delta_complex, k_I, k_f, Lambda, D, V, k_d, k_c, mu]...
                = initialize_parameters()
r_crRNA = 10^(-3);
Delta_crRNA = 10^(-3);
r_Cas9 = 10^(-3);
Delta_Cas9 = 10^(-3);
Delta_complex = 10^(-3);
k_f = .008; % coefficient of formation 
k_I = .1;
Lambda = 4;
D = 4;
V = 1;
mu = .000125;
k_d = .65;
k_c = 2;
end