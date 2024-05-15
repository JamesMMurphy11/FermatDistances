% Simple script to plot the geodesic sphere for Fermat distance.

%clear, clc, close all
x = [0;0]; % base point
num_geodesics = 100;
p=3;
rad=0.35;

xgrid = -0.5:.01:0.5;
ygrid = -0.5:.01:0.5;
density_grid = zeros(length(xgrid),length(ygrid));
for j=1:length(xgrid)
    for i=1:length(ygrid)
        density_grid(i,j) = 1 + xgrid(j);
    end
end
hf = figure;
imagesc(xgrid, ygrid, density_grid)
colorbar
alpha(0.5)

hold on
X=zeros(100,100);
Y=zeros(100,100);
for n = 1:num_geodesics
    theta = 2*pi*n/num_geodesics;
    b = [cos(theta);sin(theta)];
    y0 = [x;b]; % initial condition for ODE solver.
    [t,y] = ode45(@geodesic, [0 rad], y0);
    X(n,1:length(t))=y(:,1);
    Y(n,1:length(t))=y(:,2);
    plot(y(:,1),y(:,2))
end

xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
%title(['Fermat Geodesics of length ',num2str(rad),' for $p=$ ',num2str(p),' '], 'fontsize',18,'Interpreter','latex') 

%Add arrow
a = annotation('arrow');
a.Parent = hf.CurrentAxes;  % associate annotation with current axes
% now you can use data units
a.X = [0 X(1,length(t))];
a.Y = [0 0];
a.Color = [1 0 0];
a.LineWidth = 2;
a.HeadStyle = 'vback3';

%dim = [0.7 .25 0.3 0.3]; %for small rad
%dim = [0.7 .3 0.3 0.3]; %for large rad
if rad > .3
    dim = [.76 .5 0.05 0.05]; %for large rad
else
    dim = [.58 .5 0.05 0.05]; %for large rad
end
str = {'$\nabla \rho$'};
b = annotation('textbox',dim,'Interpreter','latex','String',str,'FitBoxToText','on','Fontsize',20,'EdgeColor','none');

