% JLS, 4/30/2018
% fplot_particles.m
% A function to plot the particles on an LxL lattice

function fplot_particles(L,Np,xpos,ypos,Nfig)

figure(Nfig); clf
plot(xpos,ypos,'ro','MarkerSize',20,'MarkerFaceColor','r')
xlabel('x'); 
ylabel('y');
axis equal
axis([0.5 L+0.5 0.5 L+0.5])
grid on
axis square
set(gca,'Ytick',1:L);
set(gca,'Xtick',1:L);
set(gca,'FontSize',14)

