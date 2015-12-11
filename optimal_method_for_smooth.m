%% Optimal method 

%% Initialize parameters
A = [4 2 -2; 2 6 0; -2 0 8]; b = [-8; -4; -2];
x0 = [50;50;50];
n = 3; %dimension of problem
N = 100; % number of iterations

x = zeros(n, N+1);
y = zeros(n, N+1);
x(:, 1) = x0;
y(:, 1) = x(:, 1);
Beta = norm(A);
Alpha = min(eig(A));
Kappa = Beta/Alpha;
t_prev = 0.5;

%% Iterations
for i = 1:100
    
    %% Step and project   
    x(:, i+1) = y(:, i) - (1/Beta)*calc_subgradient(A, b, y(:, i));
    temp = x(:, i+1); temp(temp<0) = 0; x(:, i+1) = temp;
    
    %% calculate next step size
    t_possible = roots([1 t_prev^2 - 1/Kappa -t_prev^2]);
    t_curr = t_possible(t_possible>0);
    
    %% Calc nxt y
    y(:, i+1) = x(:, i+1) + (t_prev*(1-t_prev)/(t_prev^2 + t_curr))*(x(:, i+1) - x(:, i));
    
end
