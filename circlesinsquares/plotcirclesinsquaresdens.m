mat = dlmread('packing.dat');
rad = 1; %unit circles, can assume 1nm
sig = 2*rad;

n = mat(:,1); %number of circles
L = 2+2./mat(:,2); %min edge length for sigma=2
L = L / sig; %min edge length for sigma=1

rhomax = n./(L.^2); %circles per L^2 or nm^2 in our case

rhomax5 = n./(L.^2)/25; %circles per L^2 or nm^2 in our case with the assumption that sigma = 5nm

semilogy(L,rhomax,'r','LineWidth',2);
hold on
semilogy(L,rhomax5,'LineWidth',2)
legend('\sigma=1nm','\sigma=5nm','Location','NorthWest');
ylabel('Maximum Surface Density [molecules/nm^2]','FontSize',16);
set(gca,'FontSize',16);
xlabel('Square edge length [nm]','FontSize',16);
set(gcf,'color','w');

tobefit = [L rhomax]; % doesn't really have to fit this
% one can easily observe that rhomax =~ 1/sig^2
