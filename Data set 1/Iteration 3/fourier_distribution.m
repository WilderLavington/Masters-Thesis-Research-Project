function eval = fourier_distribution(x, params, number_of_params, L)
    
    eval = zeros(size(x));
    eval = eval + params(1);
    params = params(2:end);
    counter = 0;
    
    for n = 1:2:number_of_params-1
        counter = counter + 1;
        eval = eval + params(n).*cos(counter*x*2*pi)./L;
        eval = eval + params(n+1).*sin(counter*x*2*pi)./L;
    end
    
    % eval = eval.^2;
    % eval = abs(eval);
end 
