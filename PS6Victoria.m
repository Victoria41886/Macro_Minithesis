%% Preliliminaries
%Define parameters
p_k = 1
sigma=0.5 % the PS recommends 1, but then we divide by 0, so ich chose 0.5
alpha=0.36
beta=0.98
pau= 1/beta -1 
delta= 0.01
A=1

% Calculate the steady state capital
K_ss= ((A*alpha)/(p_k*(pau+delta)))^(1/(1-alpha))

% Construct the grid for capital, 1000 gridpoints, lower bound = 0.9*K_ss,
% upper bound= K_ss
kgrid = linspace(0.9*K_ss, K_ss, 1000)

% Create an initial guess for the value function, here we choose a zeros
% vector as initial guess
V0= zeros(1000,1);

%% Main loop

%Pick a convergence criterion
con_crit= 0.00001

% Value function matrices
ugrid= zeros(1000,1000)

for i= 1:1000
  for j=1:1000
    ugrid(i,j)= (1/(1-sigma))*(A*kgrid(i).^alpha + p_k*(1-delta)*kgrid(i) - p_k*kgrid(j)).^(1-sigma)
  end
end











