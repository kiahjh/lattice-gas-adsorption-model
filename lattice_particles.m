% JLS, April 30, 2018, rev. 11/21/2022
% lattice_particles.m
% A program to simulate Np particles on an LxL lattice
%
% The program calls the functions
% finiconf.m            to generate the initial configuration
% fplot_particles.m    to plot configurations

clear; % clear all variables
clc;

% set parameters
L = 25; % side length of lattice
h = 5; % height of lattice
Np = 313; % number of particles
J = 1.0; % absolute value of particle-particle interaction energy
kappa = 1.0;

transitiontemp = 1.05; % reduced temperature kB*T/J

iflag = 0; % illustration flag, set to zero for faster simulations

% set Monte Carlo simulation parameters
kequilib = 15000; % number of equilibration steps
kobs = 25000; % number of production steps

% run this to just see how it works
% [c, e, C] = simulate(L, h, 1, kappa, J, Np, kobs, kequilib, 0, 1)
% c
% e
% C

% -----------

% run this to get isotherms for coverage(kappa)
kappa_vals = -2:0.05:5;
coverage_vals_1 = zeros(1, length(kappa_vals));
coverage_vals_2 = zeros(1, length(kappa_vals));
coverage_vals_3 = zeros(1, length(kappa_vals));
coverage_vals_4 = zeros(1, length(kappa_vals));
coverage_vals_5 = zeros(1, length(kappa_vals));
coverage_vals_6 = zeros(1, length(kappa_vals));
coverage_vals_7 = zeros(1, length(kappa_vals));

for i = 1:length(kappa_vals)
    % averages
    fprintf('i = %i\nkappa = %.2f\n\n', i, kappa_vals(i));
    coverage_vals_1(i) = (simulate(L, h, 0.50, kappa_vals(i), J, Np, kobs, kequilib, 0, 0) + simulate(L, h, 0.50, kappa_vals(i), J, Np, kobs, kequilib, 0, 0)) / 2;
    coverage_vals_2(i) = (simulate(L, h, 0.75, kappa_vals(i), J, Np, kobs, kequilib, 0, 0) + simulate(L, h, 0.75, kappa_vals(i), J, Np, kobs, kequilib, 0, 0)) / 2;
    coverage_vals_3(i) = (simulate(L, h, 1.00, kappa_vals(i), J, Np, kobs, kequilib, 0, 0) + simulate(L, h, 1.00, kappa_vals(i), J, Np, kobs, kequilib, 0, 0)) / 2;
    coverage_vals_4(i) = (simulate(L, h, 1.25, kappa_vals(i), J, Np, kobs, kequilib, 0, 0) + simulate(L, h, 1.25, kappa_vals(i), J, Np, kobs, kequilib, 0, 0)) / 2;
    coverage_vals_5(i) = (simulate(L, h, 1.50, kappa_vals(i), J, Np, kobs, kequilib, 0, 0) + simulate(L, h, 1.50, kappa_vals(i), J, Np, kobs, kequilib, 0, 0)) / 2;
    coverage_vals_6(i) = (simulate(L, h, 1.75, kappa_vals(i), J, Np, kobs, kequilib, 0, 0) + simulate(L, h, 1.75, kappa_vals(i), J, Np, kobs, kequilib, 0, 0)) / 2;
    coverage_vals_7(i) = (simulate(L, h, 2.00, kappa_vals(i), J, Np, kobs, kequilib, 0, 0) + simulate(L, h, 2.00, kappa_vals(i), J, Np, kobs, kequilib, 0, 0)) / 2;
end

figure(1); clf
hold on
h1 = plot(kappa_vals, coverage_vals_1, '-r', 'LineWidth', 2);
h2 = plot(kappa_vals, coverage_vals_2, '-g', 'LineWidth', 2);
h3 = plot(kappa_vals, coverage_vals_3, '-b', 'LineWidth', 2);
h4 = plot(kappa_vals, coverage_vals_4, '-c', 'LineWidth', 2);
h5 = plot(kappa_vals, coverage_vals_5, '-m', 'LineWidth', 2);
h6 = plot(kappa_vals, coverage_vals_6, '-y', 'LineWidth', 2);
h7 = plot(kappa_vals, coverage_vals_7, '-k', 'LineWidth', 2);
xlabel('kappa')
ylabel('coverage')
legend([h1(1), h2(1), h3(1), h4(1), h5(1), h6(1), h7(1)], 'T_{red} = 0.50', 'T_{red} = 0.75', 'T_{red} = 1.00', 'T_{red} = 1.25', 'T_{red} = 1.50', 'T_{red} = 1.75', 'T_{red} = 2.00', 'Location', 'northeastoutside')
grid on
title({['adsorptive surface coverage as a function of kappa, J = ', num2str(J)]})

% -----------

% run this to get coverage, heat capacity, and energy plots as function of
% kappa:
% kappa_vals = -2:0.05:5;
% coverage_vals = zeros(1, length(kappa_vals));
% energy_vals = zeros(1, length(kappa_vals));
% heatcap_vals = zeros(1, length(kappa_vals));

% for i = 1:length(kappa_vals)
%     fprintf('i = %i\nkappa = %.2f\n\n', i, kappa_vals(i));
%     [coverage_vals(i), energy_vals(i), heatcap_vals(i)] = simulate(L, h, transitiontemp, kappa_vals(i), J, Np, kobs, kequilib, 0, 0);
% end

% figure(1); clf
% plot(kappa_vals, coverage_vals, '-k', 'LineWidth', 2)
% xlabel('kappa')
% ylabel('coverage')
% grid on
% title({['adsorptive surface coverage as a function of kappa, T_{red} = ', num2str(transitiontemp), ', J = ', num2str(J)]})

% figure(2); clf
% plot(kappa_vals, energy_vals, '-b', 'LineWidth', 2)
% xlabel('kappa')
% ylabel('energy')
% grid on
% title({['avg. energy per particle as a function of kappa, T_{red} = ', num2str(transitiontemp), ', J = ', num2str(J)]})

% figure(3); clf
% plot(kappa_vals, heatcap_vals, '-r', 'LineWidth', 2)
% xlabel('kappa')
% ylabel('heat capacity')
% grid on
% title({['avg. heat capacity per particle as a function of kappa, T_{red} = ', num2str(transitiontemp), ', J = ', num2str(J)]})

% ---------

% run this to get coverage, heat capacity, and energy plots as function of
% Tred:
% Tred_vals = 0.4:0.05:5.0;
% coverage_vals = zeros(1, length(Tred_vals));
% energy_vals = zeros(1, length(Tred_vals));
% heatcap_vals = zeros(1, length(Tred_vals));
%
% for i = 1:length(Tred_vals)
%     fprintf('i = %i\nT_{red} = %.2f\n', i, Tred_vals(i));
%     if Tred_vals(i) < 0.3
%         fprintf('long equilib time (10,000,000)\n\n')
%         [coverage_vals(i), energy_vals(i), heatcap_vals(i)] = simulate(L, h, Tred_vals(i), kappa, J, Np, 100000, 10000000, 0, 0);
%     else
%         fprintf('short equilib time (5,000)\n\n')
%         [coverage_vals(i), energy_vals(i), heatcap_vals(i)] = simulate(L, h, Tred_vals(i), kappa, J, Np, 100000, 5000, 0, 0);
%     end
%
% end
%
% figure(1); clf
% plot(Tred_vals, coverage_vals, '-k', 'LineWidth', 2)
% xlabel('T_{red}')
% ylabel('coverage')
% grid on
% title({['adsorptive surface coverage as a function of T_{red}, kappa = ', num2str(kappa), ', J = ', num2str(J)]})
%
% figure(2); clf
% plot(Tred_vals, energy_vals, '-b', 'LineWidth', 2)
% xlabel('T_{red}')
% ylabel('energy')
% grid on
% title({['avg. energy per particle as a function of T_{red}, kappa = ', num2str(kappa), ', J = ', num2str(J)]})
%
% figure(3); clf
% plot(Tred_vals, heatcap_vals, '-r', 'LineWidth', 2)
% xlabel('T_{red}')
% ylabel('heat capacity')
% grid on
% title({['avg. heat capacity as a function of T_{red}, kappa = ', num2str(kappa), ', J = ', num2str(J)]})
