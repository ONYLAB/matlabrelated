mat = dlmread('packing.dat');
sig = 1;
n = mat(:,1);
L = 2+2./mat(:,2);
L = L * sig;
plot(n,L,'r','LineWidth',2);
hold on
plot(n,L*5,'LineWidth',2)
legend('\sigma=1nm','\sigma=5nm','Location','NorthWest');
xlabel('Number of spherical particles','FontSize',16);
set(gca,'FontSize',16);
ylabel('Square edge length [nm]','FontSize',16);
set(gcf,'color','w');
