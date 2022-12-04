function [coverage_ratio, energy, heatcap] = simulate(L, h, Tred, kappa, J, Np, kobs, kequilib, iflag, pflag)

    % set parameters
    A = L * L; % lattice area
    epsnei = -J; % energy of a nearest-neighbor pair
    k = -kappa;

    iskip = Np; % gap (in elementary steps) between illustrations

    % check if parameters make sense
    if Np > A
        disp('Too many particles, Np > A')
        Np
        A
        return % end the program
    end

    % set Monte Carlo simulation parameters
    krun = kobs + kequilib; % total number of Monte Carlo steps
    naccpt = 0; % number of accepted moves

    % Initialize random number generator
    seed = round(sum(1000 * clock));
    % seed = 17; % keep seed fixed for debugging
    rng(seed); % seed Matlab's random number generator

    % initialize arrays
    Epp = zeros(1, krun); % initialize the energy per particle to zero
    lattice = zeros(L, L, h); % set all lattice sites to 0 (empty)
    xpos = zeros(1, Np); % array for x-coordinates
    ypos = zeros(1, Np); % array for y-coordinates
    zpos = zeros(1, Np); % array for z-coordinates
    ip = zeros(1, L); is = zeros(1, L); % initialize tables for nearest neighbors
    % identify nearest neighbors (periodic boundary conditions)
    for ii = 1:L
        ip(ii) = ii - 1;
        is(ii) = ii + 1;
    end

    ip(1) = L;
    is(L) = 1;

    % set up an initial configuration
    [lattice, xpos, ypos, zpos, Nplaced] = finiconf(L, h, Np, lattice, xpos, ypos);

    if (Nplaced < Np)
        disp('not all particles could be placed')
        Nplaced
        Np
        return % end the program
    end

    % calculate the energy of the initial configuration
    energy = 0;

    for nn = 1:Np % loop over particles
        xn = xpos(nn);
        yn = ypos(nn);
        zn = zpos(nn);
        z_energies = 0;

        if zn ~= h
            % if particle isn't at the top, there might be a particle above it
            z_energies = z_energies + lattice(xn, yn, zn + 1);
        end

        if zn ~= 1
            % if particle isn't at the bottom, there might be a particle below it
            z_energies = z_energies + lattice(xn, yn, zn - 1);
        end

        neisum = lattice(ip(xn), yn, zn) + lattice(is(xn), yn, zn) + ...
            lattice(xn, ip(yn), zn) + lattice(xn, is(yn), zn);
        z_energies = z_energies * epsnei;

        if zn == 1
            % if particle is at the bottom (the adsorption surface), add
            % adsorption energy
            z_energies = z_energies + k;
        end

        energy = energy + 0.5 * epsnei * neisum + z_energies;
    end

    % plot the initial configuration
    if (pflag)
        fplot_particles(L, h, Np, xpos, ypos, zpos, 1);
        title({['Initial configuration of N_p = ', num2str(Np), ...
            ' particles'], ['on a ', num2str(L), 'x', num2str(L), 'x', num2str(h) ...
                ' lattice, energy = ', num2str(energy)]}, ...
            'FontSize', 14)
    end

    for kk = 1:krun % loop over Monte Carlo cycles

        for jj = 1:Np % jj counts elementary steps w/in one MC cycle
            nr = ceil(rand * Np); % pick a particle at random
            xold = xpos(nr); % old position coordinates
            yold = ypos(nr);
            zold = zpos(nr);
            neisum = lattice(ip(xold), yold, zold) + lattice(is(xold), yold, zold) + ...
                lattice(xold, ip(yold), zold) + lattice(xold, is(yold), zold);

            if zold ~= 1
                neisum = neisum + lattice(xold, yold, zold - 1);
            end

            if zold ~= h
                neisum = neisum + lattice(xold, yold, zold + 1);
            end

            ebefore = epsnei * neisum; % old energy of this particle

            if zold == 1
                % add surface interaction energy if on surface
                ebefore = ebefore + k;
            end

            lattice(xold, yold, zold) = 0; % clear the old lattice site
            % pick one of 6 directions for the move
            idir = ceil(rand * 6);

            % decide if move can be accepted
            aflag = 0; % acceptance flag, set to one if move is acceptedf

            if (idir == 1) % step right
                xnew = is(xold);
                ynew = yold;
                znew = zold;
            elseif (idir == 2) % step left
                xnew = ip(xold);
                ynew = yold;
                znew = zold;
            elseif (idir == 3) % step front
                xnew = xold;
                ynew = is(yold);
                znew = zold;
            elseif (idir == 4) % step back
                xnew = xold;
                ynew = ip(yold);
                znew = zold;
            elseif (idir == 5) % step up
                xnew = xold;
                ynew = yold;
                znew = zold + 1;
            elseif (idir == 6) % step down
                xnew = xold;
                ynew = yold;
                znew = zold - 1;
            else
                disp('no such direction to move in')
                idir
                kk
                jj
                return
            end

            if znew <= h && znew >= 1
                % only runs if znew is within the bounds of the lattice (can't
                % move outside)
                if (lattice(xnew, ynew, znew) == 0) % check if lattice site is available
                    neisum = lattice(ip(xnew), ynew, znew) + lattice(is(xnew), ynew, znew) + ...
                        lattice(xnew, ip(ynew), znew) + lattice(xnew, is(ynew), znew);

                    if znew ~= 1
                        neisum = neisum + lattice(xnew, ynew, znew - 1);
                    end

                    if znew ~= h
                        neisum = neisum + lattice(xnew, ynew, znew + 1);
                    end

                    eafter = epsnei * neisum; % new energy of this particle

                    if znew == 1
                        % add surface interaction energy if on surface
                        eafter = eafter + k;
                    end

                    % compare energy
                    delE = eafter - ebefore; % energy difference

                    if (eafter <= ebefore) % move lowers the energy, accept
                        aflag = 1;
                    else % move increases energy, calculate Boltzmann factor
                        boltz = exp(-delE / Tred); % Boltzmann factor

                        if rand < boltz % accept
                            aflag = 1;
                        end

                    end

                end

            end

            if (aflag == 1) % move is accepted update values
                naccpt = naccpt + 1;
                lattice(xnew, ynew, znew) = 1; % occupy new lattice site
                xpos(nr) = xnew; % update position arrays
                ypos(nr) = ynew;
                zpos(nr) = znew;
                energy = energy + delE; % update energy
            else % move was rejected, restore state
                lattice(xold, yold, zold) = 1; % reoccupy old lattice site
            end

            if (iflag == 1 && mod(jj, iskip) == 0)
                fplot_particles(L, h, Np, xpos, ypos, zpos, 3);
                pause(0.1)
            end

        end % closes loop over elementary steps

        Epp(kk) = energy / Np; % store the energy value for plotting

    end % closes loop over MC cycles

    % calculate the average energy per particle
    kused = kequilib + 1:krun;
    Eppavg = mean(Epp(kused));
    energy = Eppavg;
    % calculate the heat capacity
    heatcap = Np * (mean(Epp(kused) .* Epp(kused)) - Eppavg^2) / Tred^2;

    % recalculate the particle number from the lattice occupation
    Npcheck = sum(sum(sum(lattice)));

    if (Npcheck ~= Np)
        disp('final particle numbers do not agree')
        Np
        Npcheck
    end

    % recalculate the energy of the final configuration
    efin = 0;

    for nn = 1:Np % loop over particles
        xn = xpos(nn);
        yn = ypos(nn);
        zn = zpos(nn);
        neisum = lattice(ip(xn), yn, zn) + lattice(is(xn), yn, zn) + ...
            lattice(xn, ip(yn), zn) + lattice(xn, is(yn), zn);

        if zn ~= 1
            neisum = neisum + lattice(xn, yn, zn - 1);
        end

        if zn ~= h
            neisum = neisum + lattice(xn, yn, zn + 1);
        end

        efin = efin + 0.5 * epsnei * neisum;

        if zn == 1
            efin = efin + k;
        end

    end

    %     if (efin ~= energy)
    %         disp('final energy values do not agree')
    %         energy
    %         efin
    %     end

    particles_on_surface = sum(sum(lattice(:, :, 1)));

    coverage_ratio = particles_on_surface / L^2;

    if (pflag)
        % plot the final configuration
        fplot_particles(L, h, Np, xpos, ypos, zpos, 3);
        title({['Final configuration after ', num2str(krun), ' MC cycles'], ...
            ['N_p = ', num2str(Np), ', L = ', num2str(L), ', h = ', num2str(h) ...
                ', T_{red} = ', num2str(Tred), ', seed = ', num2str(seed)], ...
            ['Acceptance rate = ', num2str(100 * naccpt / (Np * krun), 3), ...
                ' %, energy = ', num2str(energy)]}, ...
                'FontSize', 14)

        % plot the energy as a function of time (in MC cycles)
            figure(4); clf
            plot(1:krun, Epp, 'b-', 'LineWidth', 1)
        %         ylim([-2 0])
            grid on
            xlabel('MC cycles')
            ylabel('Energy per particle, E/N_p (J)')
        title({['Energy trajectory, T_{red} = ', num2str(Tred), ...
                ', <E>/N_p = ', num2str(Eppavg, 3)], ...
            ['N_p = ', num2str(Np), ', L = ', num2str(L), ', seed = ', ...
                num2str(seed)]})
            set(gca, 'FontSize', 14)
            end
