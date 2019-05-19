function MRS
h = 0.5;        % integration step
n = 4 / h;      % crossing number

% discretization
x = linspace(0, 4, n + 1);

A = zeros(length(x));
b = zeros(length(x), 1);


A(1, 1) = 1;
A(end, end) = 1;
b(1) = -2;
b(end) = 30.0554;
% Three-point estimation
for i = 2:length(x) - 1
A(i, i) = (- 2 / (h ^ 2))-2;
A(i, i - 1) = 1 / h ^ 2;
A(i, i + 1) = 1 / h ^ 2;
b(i) =2*exp(x(i))*sin(x(i))+6*exp(x(i))*cos(x(i));
end
% Solve and print solution
w=zeros(length(b),1);
[yjac ]=jacobi(A,b,w);
[ylieb ]= liebman(A,b,w);
figure;
plot(x,yjac,'go',x,ylieb,'rx',x,fun(x),'b');
grid on;
title(['Finite difference method for h=',num2str(h)]);
legend('Jacobi + O(h^2)','Liebman + O(h^2)','Analytical solution', 'location','best');
[EavgJ EmaxJ]= err(x,yjac)
[EavgGS EmaxGS]= err(x,ylieb)

function [Eavg Emax]=err(x,y)
E = abs(fun(x)'-y); 
Eavg = 0;
for k = 1:length(x) - 1
Eavg = Eavg + (E(k) + E(k + 1)) / 2;
end
Eavg = Eavg / (length(x) - 1);          % average error
Emax = max(E);                          % maximum error
end
end

%% Jacobi method implementation
function  [y] = jacobi(A,b,y)
R = [];
D = diag(diag(A));
L = tril(A,-1); % Extract lower triangular part.
U = triu(A,1);  % Extract upper triangular part.
for i=1:60
y = -inv(D)*(L+U) * y + D\b;
%R(i) = norm(A*y - b)/norm(y);
end
end

%% Gauss–Seidel method (Liebmann method or the method of successive displacement)
function  [y] = liebman(A,b,y)
R = [];
D = diag(diag(A));
L = tril(A,-1);
U = triu(A,1);
for i=1:60
y = -inv(L+D)*(U) * y + (L+D)\b;
%R(i) = norm(A*y - b)/norm(y);
end
end

function d2y = fun(x)
d2y = exp(x).*sin(x)-2*exp(x).*cos(x);      % analytic function
end
