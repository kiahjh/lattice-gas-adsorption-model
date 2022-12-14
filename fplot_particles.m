% JLS, 4/30/2018
% fplot_particles.m
% A function to plot the particles on an LxL lattice

function fplot_particles(L, h, Np, xpos, ypos, zpos, Nfig)

    figure(Nfig); clf
    plot3(xpos, ypos, zpos, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal
    axis([0.5 L + 0.5 0.5 L + 0.5])
    grid on
    axis equal
    set(gca, 'Projection','perspective');
    set(gca, 'Ytick', 1:L);
    set(gca, 'Xtick', 1:L);
    set(gca, 'Ztick', 1:h);
    set(gca, 'FontSize', 14)

end
