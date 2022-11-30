% JLS, April 30, 2018, rev. 11/21/2022
% lattice_particles.m
% A program to simulate Np particles on an LxL lattice
%
% The program calls the function
% finiconf.m       to generate the initial configuration
% fplot_particles.m    to plot configurations

clear; % clear all variables

% set parameters
L = 10; % side length of lattice
A = L*L; % lattice area 
Np = 25; % number of particles
J = 1.0; % absolute value of particle-particle interaction energy
epsnei = -J; % energy of a nearest-neighbor pair
Tred = 2; % reduced temperature kB*T/J

iflag = 1; % illustration flag, set to zero for faster simulations
iskip = Np; % gap (in elementary steps) between illustrations

% check if parameters make sense
if Np > A
    disp('Too many particles, Np > A')
    Np
    A
    return  % end the program
end

% set Monte Carlo simulation parameters
kequilib = 5; % number of equilibration steps
kobs = 50; % number of production steps
krun = kobs + kequilib; % total number of Monte Carlo steps 
naccpt = 0; % number of accepted moves


% Initialize random number generator
seed = round(sum(1000*clock));
% seed = 17; % keep seed fixed for debugging
rng(seed); % seed Matlab's random number generator

% initialize arrays
Epp = zeros(1,krun); % initialize the energy per particle to zero
lattice = zeros(L,L); % set all lattice sites to 0 (empty)
xpos = zeros(1,Np); % array for x-coordinates
ypos = zeros(1,Np); % array for y-coordinates
ip = zeros(1,L); is = zeros(1,L); % initialize tables for nearest neighbors
% identify nearest neighbors (periodic boundary conditions)
for ii = 1:L
   ip(ii) = ii - 1;
   is(ii) = ii + 1;
end
ip(1) = L;
is(L) = 1;


% set up an initial configuration
[lattice,xpos,ypos,Nplaced] = finiconf(L,Np,lattice,xpos,ypos); 
if (Nplaced < Np)
    disp('not all particles could be placed')
    Nplaced
    Np
    return % end the program 
end

% calculate the energy of the initial configuration
energy = 0;
for nn = 1:Np  % loop over particles
    xn = xpos(nn);
    yn = ypos(nn);
    neisum = lattice(ip(xn),yn) + lattice(is(xn),yn) + ...
        lattice(xn,ip(yn)) + lattice(xn,is(yn)); 
    energy = energy + 0.5*epsnei*neisum;
end

% plot the initial configuration
fplot_particles(L,Np,xpos,ypos,1);
title({['Initial configuration of N_p = ',num2str(Np),...
    ' particles'],['on a ',num2str(L),'x',num2str(L),...
    ' lattice, energy = ',num2str(energy)]},...
    'FontSize',14)

for kk = 1:krun  % loop over Monte Carlo cycles
    for jj = 1:Np % jj counts elementary steps w/in one MC cycle
        nr = ceil(rand*Np); % pick a particle at random
        xold = xpos(nr); % old position coordinates
        yold = ypos(nr);
        neisum = lattice(ip(xold),yold) + lattice(is(xold),yold) + ...
            lattice(xold,ip(yold)) + lattice(xold,is(yold));
        ebefore = epsnei*neisum; % old energy of this particle 
        lattice(xold,yold) = 0; % clear the old lattice site
        % pick one of 4 directions for the move
        idir = ceil(rand*4); 
        if (idir == 1) % step right
            xnew = is(xold);
            ynew = yold;
        elseif (idir == 2) % step left
            xnew = ip(xold);
            ynew = yold;
        elseif (idir == 3) % step front
            xnew = xold;
            ynew = is(yold);
        elseif (idir == 4) % step back
            xnew = xold;
            ynew = ip(yold);
        else
            disp('no such direction to move in')
            idir
            kk
            jj
            return
        end
        % decide if move can be accepted
        aflag = 0; % acceptance flag, set to one if move is accepted
        if (lattice(xnew,ynew) == 0) % check if lattice site is available
            neisum = lattice(ip(xnew),ynew) + lattice(is(xnew),ynew) + ...
                lattice(xnew,ip(ynew)) + lattice(xnew,is(ynew));
            eafter = epsnei*neisum; % new energy of this particle
            % compare energy
            delE = eafter - ebefore; % energy difference
            if (eafter <= ebefore) % move lowers the energy, accept
                aflag = 1;
            else % move increases energy, calculate Boltzmann factor
                boltz = exp(-delE/Tred); % Boltzmann factor
                if rand < boltz  % accept
                    aflag = 1;
                end
            end
        end
        if (aflag == 1) % move is accepted update values
            naccpt = naccpt + 1;
            lattice(xnew,ynew) = 1; % occupy new lattice site
            xpos(nr) = xnew; % update position arrays
            ypos(nr) = ynew;
            energy = energy + delE; % update energy
        else % move was rejected, restore state
            lattice(xold,yold) = 1; % reoccupy old lattice site
        end
        
        if (iflag == 1 && mod(jj,iskip) == 0)
            fplot_particles(L,Np,xpos,ypos,3);
            pause(0.1)
        end
        
    end % closes loop over elementary steps
    
    Epp(kk) = energy/Np; % store the energy value for plotting
            
end % closes loop over MC cycles

% calculate the average energy per particle
Eppavg = mean(Epp);
           
% recalculate the particle number from the lattice occupation
Npcheck = sum(sum(lattice));
if (Npcheck ~= Np)
    disp('final particle numbers do not agree')
    Np
    Npcheck
end

% recalculate the energy of the final configuration
efin = 0;
for nn = 1:Np  % loop over particles
    xn = xpos(nn);
    yn = ypos(nn);
    neisum = lattice(ip(xn),yn) + lattice(is(xn),yn) + ...
        lattice(xn,ip(yn)) + lattice(xn,is(yn)); 
    efin = efin + 0.5*epsnei*neisum;
end
if (efin ~= energy)
    disp('final energy values do not agree')
    energy
    efin
end


% plot the final configuration
fplot_particles(L,Np,xpos,ypos,3);
title({['Final configuration after ',num2str(krun),' MC cycles'],...
    ['N_p = ',num2str(Np),', L = ',num2str(L),...
    ', T_{red} = ',num2str(Tred),', seed = ',num2str(seed)],...
    ['Acceptance rate = ',num2str(100*naccpt/(Np*krun),3),...
    ' %, energy = ',num2str(energy)]},...
    'FontSize',14)

% plot the energy as a function of time (in MC cycles)
figure(4); clf
plot(1:krun,Epp,'b-','LineWidth',1)
ylim([-2 0])
grid on
xlabel('MC cycles')
ylabel('Energy per particle, E/N_p (J)')
title({['Energy trajectory, T_{red} = ',num2str(Tred),...
    ', <E>/N_p = ',num2str(Eppavg,3)],...
    ['N_p = ',num2str(Np),', L = ',num2str(L),', seed = ',...
    num2str(seed)]})
set(gca,'FontSize',14)
