clc, clear all;
#Model Parameters
N_v = 1000; alpha_v = 0.09; alpha_c = 0.00045; Lambda_s = 20; beta_h = 0.0025;...
    mu_s = 0.02; mu_v = 0.15; mu_c = 0.02; mu_h = 1.5; alpha_h = 0.9; Lambda_H = 20.75/365/1000;...
    beta_H = 0.03; alpha_H = 0.24; mu_H = 0.000038; V_0 = 2000000; delta_H = 0.008/365;

n = 4000;
dt = zeros(1,n);
M = zeros(1,n);
N = zeros(1,n);
O = zeros(1,n);
#z = zeros(1,7);

M_0 = zeros(1,n);
N_0 = zeros(1,n);
O_0 = zeros(1,n);

M_1 = zeros(1,n);
N_1 = zeros(1,n);
O_1 = zeros(1,n);

M_2 = zeros(1,n);
N_2 = zeros(1,n);
O_2 = zeros(1,n);

M_3 = zeros(1,n);
N_3 = zeros(1,n);
O_3 = zeros(1,n);

M_4 = zeros(1,n);
N_4 = zeros(1,n);
O_4 = zeros(1,n);

M_5 = zeros(1,n);
N_5 = zeros(1,n);
O_5 = zeros(1,n);

M_6 = zeros(1,n);
N_6 = zeros(1,n);
O_6 = zeros(1,n);

#gamma=0.35; epsilon=0; kappa=0; sigma=0; v=0; eta=0; Phi =0;


h = 1;

for i = 1:n
	dt(1) = 0;
	dt(i+1) = dt(i)+h;
	M_0(1) = 13000; N_0(1) = 4000; O_0(1) = 4*10^(6);
	%Within_host model
	gamma=0; epsilon=0; kappa=0; sigma=0; v=0; eta=0; Phi =0;
	R_0w = ((1-kappa)*N_v*alpha_v*(1-epsilon)*alpha_c*Lambda_s*(1-gamma)*beta_h)/(mu_s*(alpha_v + mu_v)*((1-epsilon)*alpha_c + mu_c)*(alpha_h + mu_h));
	
	N_h = (mu_s/beta_h*(1-gamma))*(R_0w -1);
	phi = (1-exp(-(mu_H+delta_H*(1-eta))*h))/(mu_H+delta_H*(1-eta));
	
	#Modified
	M_0(i+1) = ((M_0(i) + phi*Lambda_H)*(V_0*(1+v) + O_0(i)))/((1+phi*mu_H)*(V_0*(1+v) + O_0(i)) + phi*beta_H*(1-sigma)*O_0(i));
	A = N_0(i)*((1 + phi*mu_H)*(V_0*(1+v) + O_0(i)) + beta_H*(1-sigma)*phi*O_0(i)) + (M_0(i) + phi*Lambda_H)*beta_H*(1-sigma)*phi*O_0(i);
	B = (1 + phi*mu_H + phi*delta_H*(1-eta))*((1 + phi*mu_H)*(V_0*(1+v) + O_0(i)) + beta_H*(1-sigma)*phi*O_0(i));
	N_0(i+1) = (A/B);
	X = O_0(i) + phi*N_h*alpha_h*N_0(i);
	Y = 1 + phi*alpha_H*(1+Phi);
	O_0(i+1) = (X/Y);
end


for i = 1:n
	dt(1) = 0;
	dt(i+1) = dt(i)+h;
	M(1) = 13000; N(1) = 4000; O(1) = 4*10^(6);
	%Within_host model
	gamma=0.35; epsilon=0.35; kappa=0.35; sigma=0.35; v=0.35; eta=0.175; Phi =0.175;
	R_0w = ((1-kappa)*N_v*alpha_v*(1-epsilon)*alpha_c*Lambda_s*(1-gamma)*beta_h)/(mu_s*(alpha_v + mu_v)*((1-epsilon)*alpha_c + mu_c)*(alpha_h + mu_h));
	
	N_h = (mu_s/beta_h)*(R_0w -1);
	phi = (1-exp(-(mu_H+delta_H)*h))/(mu_H+delta_H);
	
	#Modified
	M(i+1) = ((M(i) + phi*Lambda_H)*(V_0*(1+v) + O(i)))/((1+phi*mu_H)*(V_0*(1+v) + O(i)) + phi*beta_H*(1-sigma)*O(i));
	A = N(i)*((1 + phi*mu_H)*(V_0*(1+v) + O(i)) + beta_H*(1-sigma)*phi*O(i)) + (M(i) + phi*Lambda_H)*beta_H*(1-sigma)*phi*O(i);
	B = (1 + phi*mu_H + phi*delta_H*(1-eta))*((1 + phi*mu_H)*(V_0*(1+v) + O(i)) + beta_H*(1-sigma)*phi*O(i));
	N(i+1) = (A/B);
	X = O(i) + phi*N_h*alpha_h*N(i);
	Y = 1 + phi*alpha_H*(1+Phi);
	O(i+1) = (X/Y);
end

for i = 1:n
	dt(1) = 0;
	dt(i+1) = dt(i)+h;
	M_1(1) = 13000; N_1(1) = 4000; O_1(1) = 4*10^(6);
	%Within_host model
	gamma=0; epsilon=0.35; kappa=0; sigma=0; v=0; eta=0; Phi =0;
	R_0w = ((1-kappa)*N_v*alpha_v*(1-epsilon)*alpha_c*Lambda_s*(1-gamma)*beta_h)/(mu_s*(alpha_v + mu_v)*((1-epsilon)*alpha_c + mu_c)*(alpha_h + mu_h));
	N_h = (mu_s/beta_h*(1-gamma))*(R_0w -1);
	phi = (1-exp(-(mu_H+delta_H*(1-eta))*h))/(mu_H+delta_H*(1-eta));
	
	#Modified
	M_1(i+1) = ((M_1(i) + phi*Lambda_H)*(V_0*(1+v) + O_1(i)))/((1+phi*mu_H)*(V_0*(1+v) + O_1(i)) + phi*(1-sigma)*beta_H*O_1(i));
	A = N_1(i)*((1 + phi*mu_H)*(V_0*(1+v) + O_1(i)) + beta_H*(1-sigma)*phi*O_1(i)) + (M_1(i) + phi*Lambda_H)*beta_H*(1-sigma)*phi*O_1(i);
	B = (1 + phi*mu_H + phi*delta_H*(1-eta))*((1 + phi*mu_H)*(V_0*(1+v) + O_1(i)) + beta_H*(1-sigma)*phi*O_1(i));
	N_1(i+1) = (A/B);
	X = O_1(i) + phi*N_h*alpha_h*N_1(i);
	Y = 1 + phi*alpha_H*(1-Phi);
	O_1(i+1) = (X/Y);
end

for i = 1:n
	dt(1) = 0;
	dt(i+1) = dt(i)+h;
	M_2(1) = 13000; N_2(1) = 4000; O_2(1) = 4*10^(6);
	%Within_host model
	gamma=0; epsilon=0; kappa=0.35; sigma=0; v=0; eta=0; Phi =0;
	R_0w = ((1-kappa)*N_v*alpha_v*(1-epsilon)*alpha_c*Lambda_s*(1-gamma)*beta_h)/(mu_s*(alpha_v + mu_v)*((1-epsilon)*alpha_c + mu_c)*(alpha_h + mu_h));
	N_h = (mu_s/beta_h*(1-gamma))*(R_0w -1);
	phi = (1-exp(-(mu_H+delta_H*(1-eta))*h))/(mu_H+delta_H*(1-eta));
	
	#Modified
	M_2(i+1) = ((M_2(i) + phi*Lambda_H)*(V_0*(1+v) + O_2(i)))/((1+phi*mu_H)*(V_0*(1+v) + O_2(i)) + phi*(1-sigma)*beta_H*O_2(i));
	A = N_2(i)*((1 + phi*mu_H)*(V_0*(1+v) + O_2(i)) + beta_H*(1-sigma)*phi*O_2(i)) + (M_2(i) + phi*Lambda_H)*beta_H*(1-sigma)*phi*O_2(i);
	B = (1 + phi*mu_H + phi*delta_H*(1-eta))*((1 + phi*mu_H)*(V_0*(1+v) + O_2(i)) + beta_H*(1-sigma)*phi*O_2(i));
	N_2(i+1) = (A/B);
	X = O_2(i) + phi*N_h*alpha_h*N_2(i);
	Y = 1 + phi*alpha_H*(1-Phi);
	O_2(i+1) = (X/Y);
end

for i = 1:n
	dt(1) = 0;
	dt(i+1) = dt(i)+h;
	M_3(1) = 13000; N_3(1) = 4000; O_3(1) = 4*10^(6);
	%Within_host model
	gamma=0; epsilon=0; kappa=0; sigma=0.35; v=0; eta=0; Phi =0;
	R_0w = ((1-kappa)*N_v*alpha_v*(1-epsilon)*alpha_c*Lambda_s*(1-gamma)*beta_h)/(mu_s*(alpha_v + mu_v)*((1-epsilon)*alpha_c + mu_c)*(alpha_h + mu_h));
	N_h = (mu_s/beta_h*(1-gamma))*(R_0w -1);
	phi = (1-exp(-(mu_H+delta_H*(1-eta))*h))/(mu_H+delta_H*(1-eta));
	
	#Modified
	M_3(i+1) = ((M_3(i) + phi*Lambda_H)*(V_0*(1+v) + O_3(i)))/((1+phi*mu_H)*(V_0*(1+v) + O_3(i)) + phi*(1-sigma)*beta_H*O_3(i));
	A = N_3(i)*((1 + phi*mu_H)*(V_0*(1+v) + O_3(i)) + beta_H*(1-sigma)*phi*O_3(i)) + (M_3(i) + phi*Lambda_H)*beta_H*(1-sigma)*phi*O_3(i);
	B = (1 + phi*mu_H + phi*delta_H*(1-eta))*((1 + phi*mu_H)*(V_0*(1+v) + O_3(i)) + beta_H*(1-sigma)*phi*O_3(i));
	N_3(i+1) = (A/B);
	X = O_3(i) + phi*N_h*alpha_h*N_3(i);
	Y = 1 + phi*alpha_H*(1-Phi);
	O_3(i+1) = (X/Y);
end

for i = 1:n
	dt(1) = 0;
	dt(i+1) = dt(i)+h;
	M_4(1) = 13000; N_4(1) = 4000; O_4(1) = 4*10^(6);
	%Within_host model
	gamma=0; epsilon=0; kappa=0; sigma=0; v=0.35; eta=0; Phi =0;
	R_0w = ((1-kappa)*N_v*alpha_v*(1-epsilon)*alpha_c*Lambda_s*(1-gamma)*beta_h)/(mu_s*(alpha_v + mu_v)*((1-epsilon)*alpha_c + mu_c)*(alpha_h + mu_h));
	N_h = (mu_s/beta_h*(1-gamma))*(R_0w -1);
	phi = (1-exp(-(mu_H+delta_H*(1-eta))*h))/(mu_H+delta_H*(1-eta));
	
	#Modified
	M_4(i+1) = ((M_4(i) + phi*Lambda_H)*(V_0*(1+v) + O_4(i)))/((1+phi*mu_H)*(V_0*(1+v) + O_4(i)) + phi*(1-sigma)*beta_H*O_4(i));
	A = N_4(i)*((1 + phi*mu_H)*(V_0*(1+v) + O_4(i)) + beta_H*(1-sigma)*phi*O_4(i)) + (M_4(i) + phi*Lambda_H)*beta_H*(1-sigma)*phi*O_4(i);
	B = (1 + phi*mu_H + phi*delta_H*(1-eta))*((1 + phi*mu_H)*(V_0*(1+v) + O_4(i)) + beta_H*(1-sigma)*phi*O_4(i));
	N_4(i+1) = (A/B);
	X = O_4(i) + phi*N_h*alpha_h*N_4(i);
	Y = 1 + phi*alpha_H*(1-Phi);
	O_4(i+1) = (X/Y);
end

for i = 1:n
	dt(1) = 0;
	dt(i+1) = dt(i)+h;
	M_5(1) = 13000; N_5(1) = 4000; O_5(1) = 4*10^(6);
	%Within_host model
	gamma=0; epsilon=0; kappa=0; sigma=0; v=0; eta=0.175; Phi =0;
	R_0w = ((1-kappa)*N_v*alpha_v*(1-epsilon)*alpha_c*Lambda_s*(1-gamma)*beta_h)/(mu_s*(alpha_v + mu_v)*((1-epsilon)*alpha_c + mu_c)*(alpha_h + mu_h));
	N_h = (mu_s/beta_h*(1-gamma))*(R_0w -1);
	phi = (1-exp(-(mu_H+delta_H*(1-eta))*h))/(mu_H+delta_H*(1-eta));
	
	#Modified
	M_5(i+1) = ((M_5(i) + phi*Lambda_H)*(V_0*(1+v) + O_5(i)))/((1+phi*mu_H)*(V_0*(1+v) + O_5(i)) + phi*(1-sigma)*beta_H*O_5(i));
	A = N_5(i)*((1 + phi*mu_H)*(V_0*(1+v) + O_5(i)) + beta_H*(1-sigma)*phi*O_5(i)) + (M_5(i) + phi*Lambda_H)*beta_H*(1-sigma)*phi*O_5(i);
	B = (1 + phi*mu_H + phi*delta_H*(1-eta))*((1 + phi*mu_H)*(V_0*(1+v) + O_5(i)) + beta_H*(1-sigma)*phi*O_5(i));
	N_5(i+1) = (A/B);
	X = O_5(i) + phi*N_h*alpha_h*N_5(i);
	Y = 1 + phi*alpha_H*(1-Phi);
	O_5(i+1) = (X/Y);
end

for i = 1:n
	dt(1) = 0;
	dt(i+1) = dt(i)+h;
	M_6(1) = 13000; N_6(1) = 4000; O_6(1) = 4*10^(6);
	%Within_host model
	gamma=0; epsilon=0; kappa=0; sigma=0; v=0; eta=0; Phi =0.175;
	R_0w = ((1-kappa)*N_v*alpha_v*(1-epsilon)*alpha_c*Lambda_s*(1-gamma)*beta_h)/(mu_s*(alpha_v + mu_v)*((1-epsilon)*alpha_c + mu_c)*(alpha_h + mu_h));
	N_h = (mu_s/beta_h*(1-gamma))*(R_0w -1);
	phi = (1-exp(-(mu_H+delta_H*(1-eta))*h))/(mu_H+delta_H*(1-eta));
	
	#Modified
	M_6(i+1) = ((M_6(i) + phi*Lambda_H)*(V_0*(1+v) + O_6(i)))/((1+phi*mu_H)*(V_0*(1+v) + O_6(i)) + phi*(1-sigma)*beta_H*O_6(i));
	A = N_6(i)*((1 + phi*mu_H)*(V_0*(1+v) + O_6(i)) + beta_H*(1-sigma)*phi*O_6(i)) + (M_6(i) + phi*Lambda_H)*beta_H*(1-sigma)*phi*O_6(i);
	B = (1 + phi*mu_H + phi*delta_H*(1-eta))*((1 + phi*mu_H)*(V_0*(1+v) + O_6(i)) + beta_H*(1-sigma)*phi*O_6(i));
	N_6(i+1) = (A/B);
	X = O_6(i) + phi*N_h*alpha_h*N_6(i);
	Y = 1 + phi*alpha_H*(1-Phi);
	O_6(i+1) = (X/Y);
end


#inter = [gamma, epsilon, kappa, sigma, v, eta, Phi];
#z = [0.35 0.35 0.35 0.35 0.35 0.35 0.35];
#G = diag(z);
#x = G(1,:)

c_1 = [249 32 2] ./ 255;

c_2 = [94 250 81] ./ 255;

c_3 = [12 195 82] ./ 255;

c_4 = [118 63 210] ./ 255;

c_5 = [1 17 181] ./ 255;

c_6 = [104 242 41] ./ 255;


figure(1)
subplot(221)
plot(dt, M_0, 'color', c_1, 'linewidth',3, dt, M, 'color', c_2, 'linewidth',3, dt, M_1, 'color', c_3, 'linewidth',3, dt, M_2, 'color', c_4, 'linewidth',3, dt, M_3, 'color', c_5, 'linewidth',3, dt, M_4, 'color', c_6, 'linewidth',3) 
xlabel('Time(days)')
ylabel('Susceptibles(Number)')
legend('baseline','\gamma =0.35','\epsilon = 0.35','\kappa= 0.35','\sigma =0.35','v =0.35')
subplot(222)
plot(dt, N_0, 'color', c_1, 'linewidth',3, dt, N, 'color', c_2, 'linewidth',3, dt, N_1, 'color', c_3, 'linewidth',3, dt, N_2, 'color', c_4, 'linewidth',3, dt, N_3, 'color', c_5, 'linewidth',3, dt, N_4, 'color', c_6, 'linewidth',3) 
xlabel('Time(days)')
ylabel('Infected(Number)')
legend('baseline','\gamma =0.35','\epsilon = 0.35','\kappa= 0.35','\sigma =0.35','v =0.35', 'location', 'SouthEast')
subplot(212)
plot(dt, O_0, 'color', c_1, 'linewidth',3, dt, O, 'color', c_2, 'linewidth',3, dt, O_1, 'color', c_3, 'linewidth',3, dt, O_2, 'color', c_4, 'linewidth',3, dt, O_3, 'color', c_5, 'linewidth',3, dt, O_4, 'color', c_6, 'linewidth',3) 
xlabel('Time(days)')
ylabel('C-Viral Load')
legend('baseline','\gamma =0.35','\epsilon = 0.35','\kappa= 0.35','\sigma =0.35','v =0.35', 'location', 'EastOutside')



