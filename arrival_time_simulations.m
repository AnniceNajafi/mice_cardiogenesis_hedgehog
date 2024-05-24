%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@author annice
%description: In the following program we solve the DE to find s(t) and
%then calculate p
%Date: 05/22/24
%Texas A&M University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load receptor values 
filename = '~/Downloads/ptch1_values.csv'; 
receptor_vals = readtable(filename);
receptor_vals = receptor_vals(1:50, :);
receptor_vals = table2array(receptor_vals);

% Params from Lorentzian fit in R
a = 0.1560915;
b = 33.47859;
c = 50.28626;
r_peak = max(receptor_vals(:, 3)');

pp=1;
y_hold(1)=1;
y_hold_sat(1)=1;

for a = max(a)*0.01:0.01:max(a)
kk = 1;
for r_peak = max(receptor_vals(:, 3)')*0.01:0.01:max(receptor_vals(:, 3)')
% Lorentzian
syms s_sym
f_sym = a * b^2 / ((s_sym - c)^2 + b^2);
dfds_sym = diff(f_sym, s_sym);


f = matlabFunction(f_sym, 'Vars', s_sym);
dfds = matlabFunction(dfds_sym, 'Vars', s_sym);


l_9_5 = 184.75 *10^-6;  %Setting to table values for now
l = @(t) 2.^(t / 2) * l_9_5;

odefun_s = @(t, s) r_peak * dfds(s);

tspan = [0 100];
s0 = 1; 
[t_s, s] = ode23(odefun_s, tspan, s0);


p = @(t) interp1(t_s, s, t, 'linear', 'extrap') ./ l(t);


composed_function_numeric = @(t) l(t) * f(p(t) * l(t));


d_composed_function_numeric = @(t) numerical_derivative(composed_function_numeric, t);

% Solve the ODE for p(t)
y0 = 0; % Initial condition for p
[t, y] = ode23(@(t, y) odefun(t, y, d_composed_function_numeric, l, r_peak), tspan, y0);

%%% Fitting logistic function to estimate saturation point
% Logistic function
logisticFun = @(b, x) b(1) ./ (1 + exp(-b(2) * (x - b(3))));
% Initial parameters [L, k, x0]
dy_dt = gradient(y, t);
[~, idx] = min(abs(y - max(y) / 2));
x0_est = t(idx);
k_est = dy_dt(idx) / (max(y) * (1 - y(idx) / max(y)));
initialParams = [max(y), k_est ,x0_est];

% Fit the model to the data
params = nlinfit(t, y, logisticFun, initialParams);
% Fitted params
L = params(1);  % Curve maximum value
k = params(2);  % Curve steepness 
x0 = params(3); % curve's midpoint / time

saturation_point = 2*x0;
%%%
y_hold(kk, pp) = max(y);
y_hold_sat(kk,pp) = saturation_point;
kk = kk+1;
end
pp = pp+1;
end

% Plots
figure (1);
%plot(t, y, 'Color', [0.48627450980392156, 0.5764705882352941, 0.7647058823529411], 'LineWidth', 4);

subplot(2, 1, 1)
%[xx, yy] = meshgrid(max(a)*0.01:0.01:max(a),max(receptor_vals(:, 3)')*0.01:0.01:max(receptor_vals(:, 3)'));

contourf(yy, xx, y_hold_sat, 10);
colormap('pink')
title('Time of Arrival');
ylabel('Parameter a - Lorentzian');
xlabel('Peak Ptch1 Expression');
colorbar;
hold on;

subplot(2, 1, 2)
contourf(yy, xx, y_hold, 10);
colormap('summer')
title('Moved Proportion');
ylabel('Parameter a - Lorentzian');
xlabel('Peak Ptch1 Expression');
hold on;
colorbar;


function dpdt = odefun(t, p, d_composed_function_numeric, l, r_peak)
    dpdt = r_peak / l(t)^3 * d_composed_function_numeric(t);
end

function dfdx = numerical_derivative(f, x)
    h = 1e-6; % Small step for numerical differentiation
    dfdx = (f(x + h) - f(x - h)) / (2 * h);
end


