% Simple Nozzle Calculator - Rao Method
%
% Written by Vinay Williams
% Started on 17/05/21

% BNCE RAO


disp("Running Rao Solver");


%Calculate length of nozzle
                           
Ln=nozzle.ratio*(sqrt(nozzle.area_ratio)-1)*nozzle.radius_throat/tan(15*pi/180);

th_1=15*pi/180;%rad, angle at first circle
th_N=40*pi/180;%rad, angle b/w throat exit and parabola

x1=-1.5*nozzle.radius_throat*sin(th_1);%x point of circle entering the nozzle
xN=0.382*nozzle.radius_throat*sin(th_N);%X point at transition to parabola


RN=-sqrt((0.382*nozzle.radius_throat)^2-xN^2)+1.382*nozzle.radius_throat;%radius at xN


tmp1=[2*RN, 1, 0; RN^2, RN, 1; nozzle.radius_exit^2, nozzle.radius_exit, 1];

tmp2=[1/tan(th_N); xN; Ln];
dd=tmp1\tmp2;

a=dd(1);
b=dd(2);
c=dd(3);


x_c1=linspace(x1,0,50); %note, in cm
x_c2=linspace(0,xN,50); %circle exiting nozzle
y3=linspace(RN,nozzle.radius_exit,100);

ufo=size(x_c1,2); ufa=size(x_c2,2); %ufe=size(x_c3,2);

y1=-((1.5*nozzle.radius_throat)^2*ones(1,ufo)-x_c1.^2).^.5+2.5*nozzle.radius_throat*ones(1,ufo);
y2=-((0.382*nozzle.radius_throat)^2*ones(1,ufa)-x_c2.^2).^.5+1.382*nozzle.radius_throat*ones(1,ufa);

x_c3=a*y3.^2+b*y3+c*ones(1,size(y3,2));

x=[x_c1,x_c2,x_c3];
y=[y1,y2,y3];

nozzle.xpoints.inner = x;
nozzle.ypoints.inner = y;

nozzle.xpoints.outer = nozzle.xpoints.inner+nozzle.wall_thickness;
nozzle.ypoints.outer = nozzle.ypoints.inner+nozzle.wall_thickness;

if plt ==1
    plot(nozzle.xpoints.inner, nozzle.ypoints.inner);
    hold on
    plot(nozzle.xpoints.inner, -nozzle.ypoints.inner);
    title("Nozzle Contour - RAO")
    ylabel("Radial Length [m]")
    xlabel("Axial Length [m]")
    legend("Wall")  
    if plt_save == 1
        saveas(gcf,'contour.png')
    end
end

clear Ln x y x_c3 y1 y2 ufo y3 x_c2 x_c1 dd a b c tmp2 tmp1 
clear RN xN x1 th_N th_1 ufa 
