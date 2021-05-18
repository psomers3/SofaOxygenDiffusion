%clear
%close all
%load('2548Zellen7448Zyklen10e-6Nachbarn.mat')
%load('C:\Users\Johanna\Nextcloud\Studium\Bachelorarbeit\FEM\2021-04-17 FEM+Matlab_updated\Diff4e-3.mat')
for c=1:length(cellarray)
    cellarray(c).pos=cellarray(c).pos*10^6;
    cellarray(c).radius=cellarray(c).radius*10^6;
end

figure
plotSphereOxygen(cellarray,5)
set(gca,'TickLabelInterpreter','latex')
xlabel('$x$ in $\mu$m','Interpreter','Latex')
ylabel('$y$ in $\mu$m','Interpreter','Latex')
zlabel('$z$ in $\mu$m','Interpreter','Latex')
%view(axes1,[218.442574754638 20.6712828615115])


%%
close all
%load('Diff4e-3.mat')
figure
plot(entwicklung(:,1),entwicklung(:,2))
hold on
%load('Diff1e-3.mat')
plot(entwicklung(:,1),entwicklung(:,2))


figure
%load('Diff4e-3.mat')
plot(entwicklung(:,1),entwicklung(:,3))
hold on
plot(entwicklung(:,1),entwicklung(:,4))
hold on
%load('Diff1e-3.mat')
plot(entwicklung(:,1),entwicklung(:,3))
hold on
plot(entwicklung(:,1),entwicklung(:,4))




figure
%load('Diff4e-3.mat')
AnteilNorm=(entwicklung(:,3)./entwicklung(:,2))
AnteilHypo=(entwicklung(:,4)./entwicklung(:,2))
plot(entwicklung(:,1),AnteilHypo)
hold on
%load('Diff1e-3.mat')
AnteilNorm=(entwicklung(:,3)./entwicklung(:,2))
AnteilHypo=(entwicklung(:,4)./entwicklung(:,2))
plot(entwicklung(:,1),AnteilHypo)