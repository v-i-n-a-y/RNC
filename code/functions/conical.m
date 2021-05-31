% Rocket Nozzle Calculator - Conical Nozzle
%
% Written by Vinay Williams
% Written on 17/02/21
%
close all
RTOD = 180/pi;
DTOR = pi/180;

y = [];
x = [];

x(1) = 0;
y(1) = nozzle.radius_throat;

y(2) = nozzle.radius_exit;
x(2) = x(1)+((y(2)-y(1))/tan(nozzle.diverging_angle*DTOR));

x1 = [];
y1 = [];

x1(1) = 0;
y1(1) = nozzle.radius_throat;

x1(2) = -((x(2)-x(1))/2);
y1(2) = y1(1)-(tan(nozzle.converging_angle*DTOR)*x1(2));

if plt == 1
    hold on
    plot(x, y, 'cyan')
    plot(x, -y, 'cyan')
    plot(x1, y1, 'r')
    plot(x1, -y1, 'r')
    h = [2,1];
    h(1) = plot(NaN,NaN,'ocyan');
    h(2) = plot(NaN,NaN,'or');
    title("Conical Nozzle")
    ylabel("Radius [m]")
    xlabel("Length [m]")
    legend(h,'diverging region','converging region');
    if plt_save ==1 
        saveas(gcf,'contour.png')
    end
end



nozzle.xpoints.outer = [];
nozzle.ypoints.outer = [];

nozzle.xpoints.outer(1) = x1(2);
nozzle.xpoints.outer(2) = x1(1);
nozzle.xpoints.outer(3) = x(2);

nozzle.ypoints.outer(1) = y1(2);
nozzle.ypoints.outer(2) = y1(1);
nozzle.ypoints.outer(3) = y(2);

nozzle.xpoints.inner = nozzle.xpoints.outer;
nozzle.ypoints.inner = nozzle.ypoints.inner + nozzle.wall_thickness;

clear x1 x y y1 
