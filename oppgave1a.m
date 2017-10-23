function u = oppgave1a
alpha = 0; beta = 1;
n = 20; % number of nodal points
x = linspace(0,1,n); % nodal points
A = zeros(n); % system matrix
b = zeros(n,1); % right-hand side
h = diff(x); % element size
for el=1:n-1 % element loop
k = el:el+1;
A(k,k) = A(k,k) + [1,-1;-1,1]/h(el);
b(k) = b(k) + h(el)/2;
end
A = A(2:end,2:end); % remove boundary conditions
b = b(2:end);
b(end) = b(end)+beta;
u = A \ b; % solve system
u = [alpha;u];
g = @(x) -1/2*x.^2+2*x;
plot(x,u)
hold on
plot(x,g(x))