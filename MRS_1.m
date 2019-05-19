function MRS
    h =0.5; % krok calkowania
    n = 4 / h; % liczba wezlow
    % Dyskretyzacja
    x = linspace(0, 4, n + 1);
    % Inicjalizacja macierzy
    A3 = zeros(length(x));
    b = zeros(length(x), 1);
    % warunki brzegowe
    A3(1, 1) = 1;
    A3(end, end) = 1;
    b(1) = -2;
    b(end) = 30.0554;
    
    % Three-point estimation
    for i = 2:length(x) - 1
        A3(i, i) = (- 2 / (h ^ 2))-2;
        A3(i, i - 1) = 1 / h ^ 2;
        A3(i, i + 1) = 1 / h ^ 2;
        b(i) =2*exp(x(i))*sin(x(i))+6*exp(x(i))*cos(x(i));
    end
    
%     % Three-point estimation for discretization of a set of points near bounder
%     A(2,2) = -2/h^2 -2;   
%     A(2,1) = 1/h^2;      
%     A(2,3) = 1/h^2;
%     b(2) = 2*exp(x(2))*sin(x(2))+6*exp(x(2))*cos(x(2));
%     A(end-1,end-1) = -2/h^2 -2;   
%     A(end-1,end) = 1/h^2;      
%     A(end-1,end-2) = 1/h^2;
%     b(end-1) = 2*exp(x(end-1))*sin(x(end-1))+6*exp(x(end-1))*cos(x(end-1));
%     
% % Five-point estimation
%     for i = 3:length(x) - 2
%         A(i, i) = (-30/(12*h^2))-2;
%         A(i, i - 1) = 16 / (12*h^2);
%         A(i, i + 1) = 16 / (12*h^2);
%         A(i, i - 2) = -1 / (12*h^2);
%         A(i, i + 2) = -1 / (12*h^2);
%         b(i) =2*exp(x(i))*sin(x(i))+6*exp(x(i))*cos(x(i));
%     end

    
    % Solve and print solution
    y = A3 \ b;
    plot(x, y, 'rx', x, fun(x), 'b');
    title(['FDM O(h^2) with h=',num2str(h)]);
    legend('Finite difference method',' Analytical solution', 'location','best');
    E = [];
    E = abs(fun(x)'-y); 
    Eavg = 0;
    for k = 1:length(x) - 1
        Eavg = Eavg + (E(k) + E(k + 1)) / 2;
    end
    Eavg = Eavg / (length(x) - 1)                       % average error
    Emax = max(E)                                       % maximum error
end
function d2y = fun(x)
    d2y = exp(x).*sin(x)-2*exp(x).*cos(x);              % analytic function
end
