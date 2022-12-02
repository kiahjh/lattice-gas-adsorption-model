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
kappa = 2.0;

Tred = 1.0; % reduced temperature kB*T/J

iflag = 0; % illustration flag, set to zero for faster simulations

% set Monte Carlo simulation parameters
kequilib = 2000000; % number of equilibration steps
kobs = 100000; % number of production steps

% [c, e, C] = simulate(L, h, 0.1, kappa, J, Np, kobs, kequilib, 0, 1)
% c
% e
% C

Tred_vals = 0.3:0.05:5.0;
coverage_vals = zeros(1, length(Tred_vals));
energy_vals = zeros(1, length(Tred_vals));
heatcap_vals = zeros(1, length(Tred_vals));

for i = 1:length(Tred_vals)
    fprintf('i = %i\nT_{red} = %.2f\n', i, Tred_vals(i));
    if Tred_vals(i) < 0.3
        fprintf('long equilib time (10,000,000)\n\n')
        [coverage_vals(i), energy_vals(i), heatcap_vals(i)] = simulate(L, h, Tred_vals(i), kappa, J, Np, kobs, 10000000, 0, 0);
    else
        fprintf('short equilib time (5,000)\n\n')
        [coverage_vals(i), energy_vals(i), heatcap_vals(i)] = simulate(L, h, Tred_vals(i), kappa, J, Np, kobs, 5000, 0, 0);
    end
    
end

figure(1); clf
plot(Tred_vals, coverage_vals, '-k', 'LineWidth', 2)
xlabel('T_{red}')
ylabel('coverage')
grid on
title({['adsorptive surface coverage as a function of T_{red}, kappa = ', num2str(kappa), ', J = ', num2str(J)]})

figure(2); clf
plot(Tred_vals, energy_vals, '-b', 'LineWidth', 2)
xlabel('T_{red}')
ylabel('energy')
grid on
title({['avg. energy per particle as a function of T_{red}, kappa = ', num2str(kappa), ', J = ', num2str(J)]})

figure(3); clf
plot(Tred_vals, heatcap_vals, '-r', 'LineWidth', 2)
xlabel('T_{red}')
ylabel('heat capacity')
grid on
title({['avg. heat capacity as a function of T_{red}, kappa = ', num2str(kappa), ', J = ', num2str(J)]})
