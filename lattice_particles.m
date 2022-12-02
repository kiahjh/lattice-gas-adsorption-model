% JLS, April 30, 2018, rev. 11/21/2022
% lattice_particles.m
% A program to simulate Np particles on an LxL lattice
%
% The program calls the function
% finiconf.m       to generate the initial configuration
% fplot_particles.m    to plot configurations

clear; % clear all variables

% set parameters
L = 25; % side length of lattice
h = 5; % height of lattice
Np = 313; % number of particles
J = 1; % absolute value of particle-particle interaction energy

Tred = 0.1; % reduced temperature kB*T/J

iflag = 0; % illustration flag, set to zero for faster simulations

% set Monte Carlo simulation parameters
kequilib = 5000; % number of equilibration steps
kobs = 50000; % number of production steps

kappa_vals = linspace(-5, 5, 30);

% just run this to see how it works
kappa = 2.0;
[coverage,heatcap] = simulate(L, h, Tred, kappa, J, Np, kobs, kequilib, 0, 1);
coverage
heatcap

% these will graph isotherms of the coverage ratio as a function of kappa

% coverage_vals_1 = zeros(length(kappa_vals));
% coverage_vals_2 = zeros(length(kappa_vals));
% coverage_vals_3 = zeros(length(kappa_vals));
% coverage_vals_4 = zeros(length(kappa_vals));
% coverage_vals_5 = zeros(length(kappa_vals));
% coverage_vals_6 = zeros(length(kappa_vals));

% for i = 1:length(kappa_vals)
%     coverage_vals_1(i) = simulate(L, h, 0, kappa_vals(i), J, Np, kobs, kequilib, 0, 0);
%     coverage_vals_2(i) = simulate(L, h, 0.25, kappa_vals(i), J, Np, kobs, kequilib, 0, 0);
%     coverage_vals_3(i) = simulate(L, h, 0.5, kappa_vals(i), J, Np, kobs, kequilib, 0, 0);
%     coverage_vals_4(i) = simulate(L, h, 0.75, kappa_vals(i), J, Np, kobs, kequilib, 0, 0);
%     coverage_vals_5(i) = simulate(L, h, 1, kappa_vals(i), J, Np, kobs, kequilib, 0, 0);
%     coverage_vals_6(i) = simulate(L, h, 2, kappa_vals(i), J, Np, kobs, kequilib, 0, 0);
% end

% figure(5); clf
% hold on
% h1 = plot(kappa_vals, coverage_vals_1, 'b-', 'LineWidth', 2)
% h2 = plot(kappa_vals, coverage_vals_2, 'g-', 'LineWidth', 2)
% h3 = plot(kappa_vals, coverage_vals_3, 'c-', 'LineWidth', 2)
% h4 = plot(kappa_vals, coverage_vals_4, 'm-', 'LineWidth', 2)
% h5 = plot(kappa_vals, coverage_vals_5, 'r-', 'LineWidth', 2)
% h6 = plot(kappa_vals, coverage_vals_5, 'k-', 'LineWidth', 2)
% legend([h1(1), h2(1), h3(1), h4(1), h5(1), h6(1)], 'T_{red} = 0', 'T_{red} = 0.25', 'T_{red} = 0.50', 'T_{red} = 0.75', 'T_{red} = 1.0', 'T_{red} = 2.0', 'Location', 'northwest')
% title({['Coverage ratio as a function of kappa for J = ', num2str(J), ' for different temperatures']})
% grid on
% xlabel('Kappa values with J=1')
% ylabel('coverage ratio of surface')
