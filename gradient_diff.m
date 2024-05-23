%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@author annice
%description: In the following program we solve the DE to find the cell
%velocity
%Date: 05/22/24
%Texas A&M University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = '~/Downloads/ptch1_values.csv'; 
receptor_vals = readtable(filename);
receptor_vals = receptor_vals(1:50, :);
receptor_vals = table2array(receptor_vals);


syms f(s) 
%Taking the estimates for parameters of Lorentzian from R
a = 0.1560915 ; b = 33.47859 ; c = 50.28626;
f(s) = a*b^2/((s-c)^2+b^2);
dfds = diff(f, s);

s_values = linspace(0, 50, 1000);

f_values = f(s_values);
dfds_values = dfds(s_values);


figure(1);
yyaxis left
ylim([0, 0.2])
plot(s_values, f_values, 'Color',[0.3137254901960784, 0.4470588235294118, 0.4823529411764706], 'LineWidth', 4);
hold on;
plot(receptor_vals(:,1), receptor_vals(:,3),'Color', [0.2862, 0.1411, 0.2431372], 'LineWidth', 4);

yyaxis right

plot(s_values, dfds_values, 'Color', [0.49019607843137253,0.0392156862745098, 0.0392156862745098], 'LineWidth', 4);
plot(receptor_vals(:,1), dfds(1:50).*receptor_vals(:,3)', 'Color', [0.48627450980392156, 0.5764705882352941, 0.7647058823529411], 'LineWidth', 4);

hold off;

xlabel('s');
ylabel('Function Value');
title('Cell Velocity');
legend('Fitted f(s)', 'Empirical r',  'Calculated df/ds', 'Calculated v(t)');
xlim([0 50])



