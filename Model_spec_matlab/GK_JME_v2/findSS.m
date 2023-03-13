function [Y K L C I varrho Pmt F Z] = findSS

Eqq1 = solve('Z = Y + betta*gam*Z','F = Y*Pmt+betta*gam*F',...
          'Z = epsilon/(epsilon-1)*F','Z, F, Pmt');
Pmt0 = Eqq1.Pmt;
delta0 = solve('delta = delta_c+b/(1+zetta)', 'delta');

C = solve('(K^alfa*L^(1-alfa))*(1-G_Y)=C+delta*K','C');

syms betta hh chi Pmt alfa Y varrho varphi

varrho = 1/(C-hh*C)-betta*hh/(C-hh*C);
L0 = exp(-log(-chi/Pmt/(-1+alfa)/Y/varrho)/(1+varphi));

%Eqq2 = solve('1/betta=Pmt*alfa*Y/K+1-delta', 'Y=K^alfa*L^(1-alfa)', 'K, Y');
%K0 = Eqq2.K;
%Y0 = Eqq2.Y;
K0 = solve('1/betta=Pmt*alfa*Y/K+1-delta', 'K');
Y0 = solve('Y=K^alfa*L^(1-alfa)', 'Y');

%%%%%%%%%%%

betta=0.99;
sig=1.0;
hh=0.815;
varphi=0.276;
zetta=7.2;
alfa=0.33;
G_over_Y=0.2;
eta_i=1.728;
epsilon=4.167;
gam=0.779;
gam_P=0.241;
kappa_pi=1.5;
kappa_y=-0.125;
rho_i=0.0;
rho_ksi=0.66;
sigma_ksi=0.05;
rho_a=0.95;
sigma_a=0.01;
rho_g=0.95;
sigma_g=0.01;
sigma_i=0.01;
kappa=10.0;
tau=0.001;
chi=3.47600242;
b=0.03510101;
delta_c=0.02071939;
G_Y = 0.2;

%%%%%%%%%%%

Pmt = eval(Pmt0);
%delta = eval(delta0);
delta = 0.025;


% Gauss-Seidel method
% initial value
L = 1.0;
K = 1.0;

for j = 1:100
  L1 = L;
  Y = eval(Y0);
  K = eval(K0);
  L = eval(L0);
  if abs(L1-L) < 0.00001
     break
  end
end
%disp(i)
%
%Y
%K
%L
C = eval(C);
I = Y*(1-G_Y) - C;
varrho = eval(varrho);
F = eval(Eqq1.F);
Z = eval(Eqq1.Z);

% X = 1/Pmt
% R = 1/betta
