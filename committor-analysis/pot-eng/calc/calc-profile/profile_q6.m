%% profile_q6

% calculate q6 profile and density profile

clear; clc; close all;

prefix = {'../frame0-339/', '../frame340-360/', '../frame361-1000/'};
nframes = [340, 20, 640];
for i = 1:numel(prefix)
    
    %  calculate
    d = dlmread([prefix{i} 'profile_q6.dat']);
    z = d(:,1);                                 % z-dist bins
    c = d(:,2);                                 % counts in each bin
    q6tot = d(:,3);                             % total value of q6 in each bin

    binwidth = z(2) - z(1);
    xyarea = 34.08*31.974;                      % xy plane area of box from frame 0
                                                % (approx same for all frames)
    nframe = nframes(i);
    rho = c*1/(2.00*nframe*binwidth*xyarea);    % density profile
                                                % factor of 2 accounts for pbc,
                                                % each side of surface
    q6 = q6tot./c;                              % average q6 value in each bin
                                                % divide total q6 by counts
    q6(isnan(q6)) = 0;

    % plot average density profile
    figure(1)
    hold on
    plot(z,rho, 'LineWidth', 2)

    % plot average q6 profile
    figure(2)
    hold on
    plot(z,q6, 'LineWidth', 2)
end

% plot average density profile
figure(1)
xlabel('$z$/\AA', 'Interpreter', 'latex')
ylabel('$\langle \rho \rangle \cdot$ \AA$^3$', 'Interpreter', 'latex')
title('$\epsilon_{SW} = 0.29$, T = 220K', 'Interpreter', 'latex')
legend({'pre-committor', 'near committor', 'post-committor'}, 'FontSize', 10)
set(gca, 'FontSize', 16)
box on

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [1.1*fig_pos(3) 1.1*fig_pos(4)];
print(fig, 'compare-densprofile', '-dpdf')

% plot average q6 profile
figure(2)
xlabel('$z$/\AA', 'Interpreter', 'latex')
ylabel('$\langle q_6 \rangle$', 'Interpreter', 'latex')
title('$\epsilon_{SW} = 0.29$, T = 220K', 'Interpreter', 'latex')
legend({'pre-committor', 'near committor', 'post-committor'}, 'FontSize', 10, 'location', 'South')
set(gca, 'FontSize', 16)
box on

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [1.1*fig_pos(3) 1.1*fig_pos(4)];
print(fig, 'compare-q6profile', '-dpdf')