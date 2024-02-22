y = sin(pi/12);
sinfun1a = @(x) 2.5*x*cos(y) - cos(5*x)*sin(y);
sinfun1b = @(x) 2.5*x*sin(y) + cos(5*x)*cos(y);

sinfun2a = @(x) x*cos(-y) - sin(x)*sin(-y);
sinfun2b = @(x) x*sin(-y) + sin(x)*cos(-y);
% sinfun2 = @(x) (1+x).*sin(x) + 2;
sinfun3 = @(x) 3*sin(-2*x) + 5.5;

xrange = linspace(0,30*pi,1000);
xrange1 = xrange(10:50)*pi/3;
xrange2 = xrange(1:120)+5*pi;

figure; hold on;
plot(sinfun1a(xrange1),sinfun1b(xrange1),'k','LineWidth',3);
plot(sinfun1a(xrange1),-sinfun1b(xrange1)+11,'k','LineWidth',3);
plot(xrange2,sinfun3(xrange2),'k','LineWidth',3);
xticks([]); ylim([-1,13])
daspect([1 1 1])