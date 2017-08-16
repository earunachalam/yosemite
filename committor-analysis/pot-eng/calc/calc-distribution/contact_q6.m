% contact_q6

% q6 distribution in the contact (1st) layer

clear; clc; close all; hold on;

prefix = {'../frame0-339/', '../frame340-360/', '../frame361-1000/'};
for i = 1:numel(prefix)
    
    % calculate
    d = dlmread([prefix{i} 'contact_q6.dat']);
    q6 = d(:,1);                                % q6 bins
    c = d(:,2);                                 % counts in each bin
    p = c./sum(c);                              % probability for each bin


    % plot
    plot(q6, p, 'LineWidth', 2)
end

xlabel('$q_6$', 'Interpreter', 'latex')
ylabel('$p$', 'Interpreter', 'latex')
title('$\epsilon_{SW} = 0.29$, T = 220K, contact layer', 'Interpreter', 'latex')
set(gca, 'FontSize', 16)
legend({'pre-committor', 'near committor', 'post-committor'}, 'FontSize', 10)
xlim([0 1])
box on

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [1.1*fig_pos(3) 1.1*fig_pos(4)];
print(fig, 'compare-contactq6', '-dpdf')