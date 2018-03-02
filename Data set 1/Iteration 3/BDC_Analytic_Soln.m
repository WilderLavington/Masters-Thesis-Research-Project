function soln = BDC_Analytic_Soln(r)
    a_0 = 0;
    N = length(r);
    if all(r)
        for j = 1:N
           prod = 1;
           for k = 1:j
               prod = prod*(r(k));
           end
           a_0 = a_0 + prod;
        end
        a_0 = a_0 + 1;
    elseif all(r)
        a_0 = 1;
    else
        a_0 = inf;
    end
    soln = 1/a_0;
end