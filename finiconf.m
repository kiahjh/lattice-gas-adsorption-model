% JLS, 4/30/2018
% finiconf.m
% A function to place Np particles on an LxL lattice

function [lattice, xpos, ypos, zpos, Nplaced] = finiconf(L, h, Np, lattice, xpos, ypos)

    nn = 0; % number of particles that have been placed = particle counter
    ntry = 0; % counter for placement attempts
    nmax = 1000; % maximum number of attempts

    while (nn < Np && ntry < nmax)
        ntry = ntry + 1; % update attempt counter
        r = rand; % create a random number between 0 and 1
        xtry = ceil(r * L); % random integer between 1 and L
        r = rand; % create a random number between 0 and 1
        ytry = ceil(r * L); % random integer between 1 and L
        r = rand; % create a random number between 0 and 1
        ztry = ceil(r * h); % random integer between 1 and h

        if (lattice(xtry, ytry, ztry) == 0) % site is empty
            nn = nn + 1; % update particle counter
            xpos(nn) = xtry; % assign x-coordinate to particle nn
            ypos(nn) = ytry; % assign y-coordinate to particle nn
            zpos(nn) = ztry;
            lattice(xtry, ytry, ztry) = 1; % declare lattice site occupied
        end

    end

    Nplaced = nn;

    if (Nplaced < Np) % could not place all particles
        return
    end

    return;

end
