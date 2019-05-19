function bvp_shooting
    close all;
    h = 0.5;                              % integration step
    x = 0:h:4;
    DESIRED_BND_VALUE = 30.0554; 			% second boundary value
    gueses = [- 0.55, 0.55] 				% state variables shots
    sol = [];
    hold on
    for g = 1:length(gueses)
        y = [-2
        gueses(g)]; % guess ith
        for i = 1:length(x) - 1
% Euler
 y(:, i + 1) = y(:, i) + h * f(x(i), y(:, i)); 
% Heun
% y(:,i+1) = y(:,i) + h/2*(f(x(i), y(:,i)) + f(x(i+1), y(:,i) + h*f(x(i),y(:,i)))); 
        end
        sol(g) = y(1, end);
        plot(x, y(1, :));
    end
    for shoot = 1:10
        g = g + 1;
        gueses(g) = spline(sol, gueses, DESIRED_BND_VALUE);
            % returns a vector of interpolated values s corresponding to the query points in                 
            % DESIRED_BND_VALUE. The values of s are determined by cubic spline interpolation 
            % of sol and queses.
        y = [-2
        gueses(g)]; % guess i-th
        for i = 1:length(x) - 1
            y(:, i + 1) = y(:, i) + h * f(x(i), y(:, i));
            % y(:,i+1) = y(:,i) + h/2*(f(x(i), y(:,i)) + f(x(i+1), y(:,i) + h * f(x(i), y(:, i))));
            E(i) = norm(y(1, i) - fun(x(i)));
        end
        sol(g) = y(1, end)
        plot(x, y(1, :));
        if (abs(sol(g) - DESIRED_BND_VALUE) < 1e-6)
             fprintf('Found solution %f for guess = %f\n', sol(g),gueses(g));
            figure
            plot(x, fun(x), '-x', x, y(1, :), '-o');
            break
        end
    end
    Emax = max(E)                                       % maximum error
    Eavg = 0;
    for i = 1:1:length(E) - 1
        Eavg = Eavg + (E(i) + E(i + 1)) / 2;
    end
    Eavg = Eavg / (length(E) - 1)                       % average error
end
function dy = f(x, y)
    dy = [y(2)
        2*exp(x)*sin(x) + 6*exp(x)*cos(x) + 2*y(1)];	% calculate the derivative
end
function d2y = fun(a)	
    d2y = exp(a).*sin(a)-2*exp(a).*cos(a); 			% analytic function
end
