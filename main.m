%% spherical harmonic expansion of gravity potential about misshapen asteroid

EquivalentDim = 140E-3; % equivalent diameter; the axial ratio of 3.07x1.0x1.0
SF = (EquivalentDim/16.84); % scaling factor
ref = 16*SF; % reference radius (km), 16km is the reference radius of Eros shape model
GM = 4.46275472004E-04*SF^3; % GM value of the target asteroid scaled down from Eros' GM value (km3/s2)
AR = [34.4 11.2 11.2]*SF; % altitude = 1/2 maximum extent+10 metres
distance =AR(1)/2+10/10^3; % (km)

potential=spharmAst(distance,0,0,ref,GM);
phi=zeros(181,361);
lam=zeros(181,361);
Potential=zeros(101,101);
for i =1:181
    phi(i,:)=-pi/2+pi*(i-1)/180;
    for j =1:361
        lam(:,j)=-pi+2*pi*(j-1)/360;
        Potential(i,j)=spharmAst(distance,phi(i,j),lam(i,j),ref,GM);
    end
end
save('dataset.mat')
%
%% Plot Gravity Potential at the given altitude
close(figure(1))
figure(1)
set(1,'Position',[40 40 1200 640],'Color',[1 1 1]);
%
load('dataset.mat')
DPH = 360/4;
plot(0:1/DPH:360/DPH,Potential(91,:)/distance*10^3)
set(gca,'FontName','Arial','FontSize',12);
xlabel('Time (hour)','FontName','Arial','FontSize',12);
ylabel('Gravity (m/s^2)','FontName','Arial','FontSize',12);
title(['Gravity change about 140m size target asteroid at minimal standoff distance=' num2str((round((distance-AR(1)/2)*1000))) 'm' ' (rotation period of target=4h)']...
    ,'FontName','Arial','FontSize',12)

%% Plot Gravity Potential at the given altitude
close(figure(2))
figure(2)
set(2,'Position',[40 40 1200 640],'Color',[1 1 1]);
%
load('dataset.mat')
xlin=linspace(-pi,pi, 25);
ylin=linspace(-pi/2,pi/2, 25);
[X,Y]=meshgrid(xlin,ylin);
Z=griddata(lam,phi,Potential,X,Y,'linear');
[C, h] = contourf(X*180/pi,Y*180/pi,Z);
set(h,'Edgecolor','none','TextStep',get(h,'LevelStep'))
clabel(C,h,'FontSize',15,'LabelSpacing',600,'Color','k','FontName','Arial')
colormap default
set(gca,'FontName','Arial','FontSize',12);
xlabel('longitude (deg)','FontName','Arial','FontSize',12);
ylabel('latitude (deg)','FontName','Arial','FontSize',12);
cbar = colorbar('FontName','Arial','FontSize',12);
ylabel(cbar,'Gravitational potential (km^2/s^2)','FontName','Arial','FontSize',12);
title(['Gravity potential map of 140m size target asteroid with minimal standoff distance=' num2str((round((distance-AR(1)/2)*1000))) 'm']...
    ,'FontName','Arial','FontSize',12)
