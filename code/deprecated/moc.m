% Rocket Nozzle Calculator - Method of Characteristics
%
% Written by Vinay Williams
% Written on 17/02/21
%
%

RTOD = 180/pi;
DTOR = pi/180;

P = []; %x axis points

% Prandtl Meyer Function
A = sqrt((nozzle.specific_heat_ratio+1)/(nozzle.specific_heat_ratio-1));
B = (nozzle.specific_heat_ratio-1)/(nozzle.specific_heat_ratio+1);
v_PM = @(x) A*atan(sqrt(B*(x^2-1))) - atan(sqrt(x^2-1));

% Diverging Angle
nozzle.moc.theta_max = 0.5*v_PM(nozzle.exit_mach)*RTOD; % Max Angle
nozzle.moc.delta_theta = (90-nozzle.moc.theta_max) - fix(90-nozzle.moc.theta_max); % 
nozzle.moc.angle_increment = nozzle.moc.delta_theta*DTOR;
n = nozzle.moc.theta_max*2;

for m = 2:n+1
    nozzle.temperature_ambient(m) = (nozzle.moc.delta_theta + (m-1))*DTOR;
    x_int = [1 1.0001*nozzle.exit_mach];
    func = @(x) nozzle.temperature_ambient(m) - v_PM(x);
    M(m) = fzero(func,x_int);
    P(m) = 0 + nozzle.radius_throat*tan(nozzle.temperature_ambient(m)); %X-AXIS POINnozzle.temperature_ambientS
    RR(m) = -nozzle.radius_throat/P(m);
    LR(m) = tan(nozzle.temperature_ambient(m)+asin(1/M(m)));
    SL(m) = -RR(m);
end

P(1) = [];
l = length(P);

    for j = 1:l
        P1 = [0 nozzle.radius_throat];
        P2 = [P(j) 0];
        
        if plt ==1 
            plot(P2,P1,'k')
            plot(P2,-P1,'k')
            hold on
            xlabel('Axial Length [m]')
            ylabel('Radial Length [m]')
        end
    end

LR(1) = []; RR(1) = [];
SL(1) = [];
F = RR(m-1);

    for c = 1:length(P)-1
        x(c) = (nozzle.radius_throat+SL(c)*P(c))/(SL(c)-F);
        y(c) = F*x(c)+nozzle.radius_throat;
        X_P = [P(c) x(c)];
        Y_P = [0 y(c)];
        if plt ==1 
            plot(X_P,Y_P,'b');
            plot(X_P,-Y_P,'b');
            hold on
        end
    end

nozzle.temperature_ambientM = nozzle.moc.theta_max*DTOR;

    xw(1) = (nozzle.radius_throat+SL(1)*P(1))/(SL(1)-tan(nozzle.temperature_ambientM));
    yw(1) = tan(nozzle.temperature_ambientM)*xw(1)+nozzle.radius_throat;
    X_P2 = [P(1) xw];
    Y_P2 = [P(2) yw];
    if plt == 1
        plot(X_P2,Y_P2,'g');
        plot(X_P2,-Y_P2,'g');
    end

%DIVIDE (delta slopes)
nozzle.moc.delta_thetaW = tan(nozzle.temperature_ambientM)/(length(P)-1);
s(1) = tan(nozzle.temperature_ambientM);
b(1) = nozzle.radius_throat;

    for k = 2:length(P)-1
        s(k) = tan(nozzle.temperature_ambientM)-(k-1)*nozzle.moc.delta_thetaW; %slope
        b(k) = yw(k-1)-s(k)*xw(k-1); %y-int
        xw(k) = (b(k)+SL(k)*P(k))/(SL(k)-s(k));
        yw(k) = s(k)*xw(k)+b(k);
        X_P3 = [x(k) xw(k)];
        Y_P3 = [y(k) yw(k)];
        if plt == 1
            plot(X_P3,Y_P3,'r');
            plot(X_P3,-Y_P3,'r');
            hold on
        end
    end
    
    
    xf = (b(length(b))+SL(length(SL))*P(length(P)))/SL(length(SL));
    yf = b(length(b));
    X_F = [P(length(P)) xf];
    Y_F = [0 yf];
%     plot(X_F,Y_F,'r');
%     plot(X_F,-Y_F,'r');
%     
    xw = [0 xw];
    yw = [nozzle.radius_throat yw];
  
    if plt ==1
    plot(xw, yw, 'cyan')
    plot(xw, -yw, 'cyan')
    h = plot(NaN,NaN,'ocyan');
    legend(h,'wall');
    title("Bell Nozzle - Method of Characteristics")
    if plt_save ==1 
        saveas(gcf,'contour.png')
    end
    end
    nozzle.xpoints.outer = xw;
    nozzle.ypoints.outer = yw+nozzle.wall_thickness;
    
    nozzle.xpoints.inner = xw;
    nozzle.ypoints.inner = yw;
    


clear A b B c DTOR RTOD f func h j k l LR m M n P P1 P2 RR s SL v_PM
clear x_int X_P2 X_P3 y Y_P2 Y_P3 yf x X_F x_w y_w Y_F Y_P yw xw F X_P xf
