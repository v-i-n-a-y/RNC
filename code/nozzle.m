classdef nozzle
    
    % Initally defined variables
    properties (Dependent = false)
        force
        impulse
        mass_flow_rate_exhaust
        y                        % Specific Heat Ratio
        chamber_pressure
        molar_mass_exhaust
        chamber_temperature
        altitude
        wall_thickness
        kind
        diverging_angle
        converging_angle
        ratio
    end
    
    properties (Constant, Access = private)
        g = 9.812;                % Gravitational Acceleration
        R = 8.314;                % Universal Gas Constant
        bar = 100000;             % Bar to pascal
        rtod = 180/pi;            % Radians to degrees
        dtor = pi/180;            % Degrees to radians
    end
    
    % Calculated Variables
    properties (Dependent = true)
        ambient_temperature
        ambient_pressure
        burn_time
        specific_impulse
        weight_flow_rate
        throat_temperature
        throat_pressure
        throat_area
        throat_radius
        exit_mach
        exit_area
        exit_radius
        exit_pressure
        exit_temperature
        exit_speed_sound
        exit_velocity
        pressure_ratio
        temperature_ratio
        pressure_ratio_2
        throat_velocity
        area_ratio
        R_dot                    % Gas Constant of Exhaust
        contour
        
    end
    
    methods (Static)
        function plot(thenozzle)
            close all
            hold on
            
            contour1 = thenozzle.contour; 
            
            plot(contour1.x.inner, contour1.y.inner);
            plot(contour1.x.inner, -contour1.y.inner);
        end
    end
    
    
    methods
        
        % Constructor
        function thisnozzle = nozzle(force, impulse, ...
                mass_flow_rate_exhaust, y, chamber_pressure, ...
                molar_mass_exhaust,chamber_temperature,altitude, ...
                wall_thickness,kind, converging_angle, diverging_angle, ratio)
            
            thisnozzle.force = force;
            thisnozzle.impulse = impulse;
            thisnozzle.chamber_temperature = chamber_temperature;
            thisnozzle.altitude = altitude * 1000;
            thisnozzle.kind = kind;
            thisnozzle.molar_mass_exhaust = molar_mass_exhaust/1000;
            thisnozzle.mass_flow_rate_exhaust = mass_flow_rate_exhaust;
            thisnozzle.wall_thickness = wall_thickness;
            thisnozzle.chamber_pressure = chamber_pressure*100000;
            thisnozzle.y = y;
            
            if kind == 3
                if ~exist('converging_angle', 'var')
                    error("Conical Nozzle chosen: Please input converging angle")
                elseif ~exist('diverging_angle', 'var')
                    error("Conical Nozzle chosen: Please input diverging angle")
                    
                else
                    thisnozzle.diverging_angle = diverging_angle;
                    thisnozzle.converging_angle = converging_angle;
                end
            elseif kind == 1
                if ~exist('converging_angle', 'var')
                    error("Rao Nozzle chosen: Please input converging angle")
                    
                elseif ~exist('diverging_angle', 'var')
                    error("Rao Nozzle chosen: Please input diverging angle")
                    
                elseif ~exist('ratio', 'var')
                    error("Rao Nozzle chosen: Please input ratio")
                else
                    thisnozzle.diverging_angle = diverging_angle;
                    thisnozzle.converging_angle = converging_angle;
                    thisnozzle.ratio = ratio; 
                end
            end
        end
        
        % CALCULATIONS
        
        % Gas Constant
        function R_dot = get.R_dot(thenozzle)
            R_dot = thenozzle.R/thenozzle.molar_mass_exhaust;
        end
        
        % Calculate ambient temperature
        function ambient_temperature = get.ambient_temperature(thenozzle)
            if (11000>thenozzle.altitude) && (thenozzle.altitude<25000)
                ambient_temperature = -56.46;
            elseif thenozzle.altitude>=25000
                ambient_temperature = -131.21 + 0.00299*thenozzle.altitude;
            else
                ambient_temperature = 15.04 - 0.00649*thenozzle.altitude;
                
            end
        end
        
        % Calculate ambient pressure
        function ambient_pressure = get.ambient_pressure(thenozzle)
            if (11000>thenozzle.altitude) && (thenozzle.altitude<25000)
                ambient_pressure = 1000*(22.65*exp(1.73-0.000157*thenozzle.altitude));
            elseif thenozzle.altitude>=25000
                ambient_pressure = 1000*(2.488*((thenozzle.ambient_temperature+273.1)/216.6)^-11.388);
            else
                ambient_pressure = 1000*(101.29*((thenozzle.ambient_temperature+273.1)/288.08)^5.256);
            end
        end
        
        
        % Burn time
        function burn_time = get.burn_time(thenozzle)
            burn_time = thenozzle.impulse/thenozzle.force;
        end
        
        % Specific Impulse
        function specific_impulse = get.specific_impulse(thenozzle)
            specific_impulse = thenozzle.force/thenozzle.mass_flow_rate_exhaust;
        end
        
        % Weight Flow Rate
        function weight_flow_rate = get.weight_flow_rate(thenozzle)
            weight_flow_rate = thenozzle.mass_flow_rate_exhaust * thenozzle.g;
        end
        
        % Throat Temperature
        function throat_temperature = get.throat_temperature(thenozzle)
            throat_temperature = thenozzle.chamber_temperature*(1/(1+((thenozzle.y-1)/2)));
        end
        
        % Throat Pressure
        function throat_pressure = get.throat_pressure(thenozzle)
            throat_pressure = thenozzle.chamber_pressure*(1+((thenozzle.y-1)/2))^(-thenozzle.y/(thenozzle.y-1));
        end
        
        % Throat Area
        function throat_area = get.throat_area(thenozzle)
            throat_area = (thenozzle.mass_flow_rate_exhaust/thenozzle.throat_pressure)*sqrt((thenozzle.R_dot*thenozzle.throat_temperature)/(thenozzle.y*thenozzle.g));
        end
        
        % Throat Radius
        function throat_radius = get.throat_radius(thenozzle)
            throat_radius = sqrt(thenozzle.throat_area/pi);
        end
        
        % Exit Mach
        function exit_mach = get.exit_mach(thenozzle)
            exit_mach = sqrt((2/(thenozzle.y-1))*((thenozzle.chamber_pressure/thenozzle.ambient_pressure)^((thenozzle.y-1)/thenozzle.y)-1));
        end
        
        % Exit Area
        function exit_area = get.exit_area(thenozzle)
            y1 = thenozzle.y-1;
            y2 = thenozzle.y+1;
            exit_area = (thenozzle.throat_area/thenozzle.exit_mach) *((1+(y1/2)*thenozzle.exit_mach^2)/(y2/2))^(y2/(2*y1));
        end
        
        % Exit Radius
        function exit_radius = get.exit_radius(thenozzle)
            exit_radius = sqrt(thenozzle.exit_area/pi);
        end
        
        % Exit Pressure
        function exit_pressure = get.exit_pressure(thenozzle)
            exit_pressure = thenozzle.throat_pressure *(1+((thenozzle.y-1)/2)* (thenozzle.exit_mach^2))^(-thenozzle.y/(thenozzle.y-1));
        end
        
        % Exit Temperature
        function exit_temperature = get.exit_temperature(thenozzle)
            exit_temperature = thenozzle.chamber_temperature/(1+((thenozzle.y-1)/2)*thenozzle.exit_mach^2);
        end
        
        % Exit Speed of Sound
        function exit_speed_sound = get.exit_speed_sound(thenozzle)
            exit_speed_sound = sqrt(thenozzle.y*thenozzle.R_dot*thenozzle.exit_temperature);
        end
        
        % Exit Velocity
        function exit_velocity = get.exit_velocity(thenozzle)
            exit_velocity = thenozzle.exit_speed_sound * thenozzle.exit_mach;
        end
        
        % Pressure Ratio
        function pressure_ratio = get.pressure_ratio(thenozzle)
            pressure_ratio = thenozzle.ambient_pressure/thenozzle.chamber_pressure;
        end
        
        % Temperature Ratio
        function temperature_ratio = get.temperature_ratio(thenozzle)
            temperature_ratio = (2*thenozzle.y*thenozzle.R_dot*thenozzle.chamber_temperature)/(thenozzle.y-1);
        end
        
        % Pressure Ratio 2
        function pressure_ratio_2 = get.pressure_ratio_2(thenozzle)
            pressure_ratio_2 = thenozzle.ambient_pressure/thenozzle.chamber_pressure^((thenozzle.y-1)/thenozzle.y);
        end
        
        % Throat Velocity
        function throat_velocity = get.throat_velocity(thenozzle)
            throat_velocity = sqrt((2*thenozzle.y*thenozzle.R_dot*thenozzle.chamber_temperature)/(thenozzle.y+1));
        end
        
        % Area Ratio
        function area_ratio = get.area_ratio(thenozzle)
            area_ratio = thenozzle.exit_area/thenozzle.throat_area;
        end
        
        % Contour
        function contour = get.contour(thenozzle)
            if thenozzle.kind == 1
                disp("CD Nozzle: Rao's Method")
                contour = rao(thenozzle);
            elseif thenozzle.kind == 2
                disp("CD Nozzle: Method of Characteristics")
                contour = moc(thenozzle);
            elseif thenozzle.kind == 3
                disp("CD Nozzle: Conical")
                contour = conical(thenozzle);
            elseif thenozzle.kind == 4
                disp("Aerospike: Radial")
            elseif thenozzle.kind == 5
                disp("Aerospike: Linear")
            else
                disp("Wrong Input")
            end
            
        end
        
        % Method of characteristics method
        function ctr = moc(thenozzle)
            
            % Prandtl Meyer Function
            A = sqrt((thenozzle.y+1)/(thenozzle.y-1));
            B = (thenozzle.y-1)/(thenozzle.y+1);
            v_PM = @(x) A*atan(sqrt(B*(x^2-1))) - atan(sqrt(x^2-1));
            
            % Diverging Angle
            theta_max = 0.5*v_PM(thenozzle.exit_mach)*thenozzle.rtod; % Max Angle
            delta_theta = (90-theta_max) - fix(90-theta_max); %
            n = theta_max*2;
            
            [ambienttemp, M, P, RR, LR, SL, x, y0, s, b, xw, yw] = deal(0);
            
            for m = 2:n+1
                
                ambienttemp(m) = (delta_theta + (m-1))*thenozzle.dtor;
                
                x_int = [1 1.0001*thenozzle.exit_mach];
                func = @(x) ambienttemp(m) - v_PM(x);
                M(m) = fzero(func,x_int);
                P(m) = 0 + thenozzle.throat_radius*tan(ambienttemp(m)); %X-AXIS POINT
                RR(m) = -thenozzle.throat_radius/P(m);
                LR(m) = tan(ambienttemp(m)+asin(1/M(m)));
                SL(m) = -RR(m);
            end
            
            P(1) = [];
            l = length(P);
            for j = 1:l
                P1 = [0 thenozzle.throat_radius];
                P2 = [P(j) 0];
                
%                 plot(P2,P1,'k')
%                 plot(P2,-P1,'k')
%                 hold on
%                 xlabel('Axial Length [m]')
%                 ylabel('Radial Length [m]')
                
            end
            
            
    
            RR(1) = [];
            SL(1) = [];
            F = RR(m-1);
            
            for c = 1:length(P)-1
                x(c) = (thenozzle.throat_radius+SL(c)*P(c))/(SL(c)-F);
                y0(c) = F*x(c)+thenozzle.throat_radius;
                X_P = [P(c) x(c)];
                Y_P = [0 y0(c)];
                
%                 plot(X_P,Y_P,'b');
%                 plot(X_P,-Y_P,'b');
%                 hold on
               
            end
            
            ambient_temperatureM = theta_max*thenozzle.dtor;
            
            xw(1) = (thenozzle.throat_radius+SL(1)*P(1))/(SL(1)-tan(ambient_temperatureM));
            yw(1) = tan(ambient_temperatureM)*xw(1)+thenozzle.throat_radius;
            X_P 2 = [P(1) xw];
            Y_P2 = [P(2) yw];
            
%             plot(X_P2,Y_P2,'g');
%             plot(X_P2,-Y_P2,'g');
            
            
            %DIVIDE (delta slopes)
            delta_thetaW = tan(ambient_temperatureM)/(length(P)-1);
            s(1) = tan(ambient_temperatureM);
            b(1) = thenozzle.throat_radius;
            
            
            for k = 2:length(P)-1
                s(k) = tan(ambient_temperatureM)-(k-1)*delta_thetaW; %slope
                b(k) = yw(k-1)-s(k)*xw(k-1); %y-int
                xw(k) = (b(k)+SL(k)*P(k))/(SL(k)-s(k));
                yw(k) = s(k)*xw(k)+b(k);
                X_P3 = [x(k) xw(k)];
                Y_P3 = [y0(k) yw(k)];
                
%                 plot(X_P3,Y_P3,'r');
%                 plot(X_P3,-Y_P3,'r');
%                 hold on
                
            end
            
            
            xf = (b(length(b))+SL(length(SL))*P(length(P)))/SL(length(SL));
            yf = b(length(b));
            X_F = [P(length(P)) xf];
            Y_F = [0 yf];
%             plot(X_F,Y_F,'r');
%             plot(X_F,-Y_F,'r');
            
            xw = [0 xw];
            yw = [thenozzle.throat_radius yw];
            
            
%             plot(xw, yw, 'cyan')
%             plot(xw, -yw, 'cyan')
%             h = plot(NaN,NaN,'ocyan');
%             legend(h,'wall');
%             title("Bell Nozzle - Method of Characteristics")
           
            ctr.x.outer = xw;
            ctr.y.outer = yw+thenozzle.wall_thickness;
            
            ctr.x.inner = xw;
            ctr.y.inner = yw;
           
        end
        
        % Rao's method method
        function ctr = rao(thenozzle)
            
            %Calculate length of nozzle
            
            Ln=thenozzle.ratio*(sqrt(thenozzle.area_ratio)-1)*thenozzle.throat_radius/tan(15*pi/180);
            
            th_1=15*pi/180;%rad, angle at first circle
            th_N=40*pi/180;%rad, angle b/w throat exit and parabola
            
            x1=-1.5*thenozzle.throat_radius*sin(th_1);%x point of circle entering the nozzle
            xN=0.382*thenozzle.throat_radius*sin(th_N);%X point at transition to parabola
            
            
            RN=-sqrt((0.382*thenozzle.throat_radius)^2-xN^2)+1.382*thenozzle.throat_radius;%radius at xN
            
            
            tmp1=[2*RN, 1, 0; RN^2, RN, 1; thenozzle.exit_radius^2, thenozzle.exit_radius, 1];
            
            tmp2=[1/tan(th_N); xN; Ln];
            dd=tmp1\tmp2;
            
            a=dd(1);
            b=dd(2);
            c=dd(3);
            
            
            x_c1=linspace(x1,0,50); %note, in cm
            x_c2=linspace(0,xN,50); %circle exiting nozzle
            y3=linspace(RN,thenozzle.exit_radius,100);
            
            ufo=size(x_c1,2); ufa=size(x_c2,2); %ufe=size(x_c3,2);
            
            y1=-((1.5*thenozzle.throat_radius)^2*ones(1,ufo)-x_c1.^2).^.5+2.5*thenozzle.throat_radius*ones(1,ufo);
            y2=-((0.382*thenozzle.throat_radius)^2*ones(1,ufa)-x_c2.^2).^.5+1.382*thenozzle.throat_radius*ones(1,ufa);
            
            x_c3=a*y3.^2+b*y3+c*ones(1,size(y3,2));
            
            x=[x_c1,x_c2,x_c3];
            y0=[y1,y2,y3];
            
            ctr.x.inner = x;
            ctr.y.inner = y0;
            
            ctr.x.outer = ctr.x.inner + thenozzle.wall_thickness;
            ctr.y.outer = ctr.y.inner + thenozzle.wall_thickness;
   
            plot(ctr.x.inner, ctr.y.inner);
            hold on
            plot(ctr.x.inner, -ctr.y.inner);
            title("Nozzle Contour - RAO")
            ylabel("Radial Length [m]")
            xlabel("Axial Length [m]")
            legend("Wall")
       
        end
        
        % Conical Nozzle method
        function ctr = conical(thenozzle)
            
            
            y0 = [];
            x = [];
  
            x(1) = 0;
            y0(1) = thenozzle.throat_radius;
            
            y0(2) = thenozzle.exit_radius;
            x(2) = x(1)+((y0(2)-y0(1))/tan(thenozzle.diverging_angle*thenozzle.dtor));
          
            x1 = [];
            y1 = [];
            
            x1(1) = 0;
            y1(1) = thenozzle.throat_radius;
            
            x1(2) = -((x(2)-x(1))/2);
            y1(2) = y1(1)-(tan(thenozzle.converging_angle*thenozzle.dtor)*x1(2));
  
            
            hold on
            plot(x, y0, 'cyan')
            plot(x, -y0, 'cyan')
            plot(x1, y1, 'r')
            plot(x1, -y1, 'r')
            h = [2,1];
            h(1) = plot(NaN,NaN,'ocyan');
            h(2) = plot(NaN,NaN,'or');
            title("Conical Nozzle")
            ylabel("Radius [m]")
            xlabel("Length [m]")
            legend(h,'diverging region','converging region');

            
            
            
            ctr.x.outer = [];
            ctr.y.outer = [];
            
            ctr.x.outer(1) = x1(2);
            ctr.x.outer(2) = x1(1);
            ctr.x.outer(3) = x(2);
            
            ctr.y.outer(1) = y1(2);
            ctr.y.outer(2) = y1(1);
            ctr.y.outer(3) = y0(2);
            
            ctr.x.inner = ctr.x.outer;
            ctr.y.inner = ctr.y.outer + thenozzle.wall_thickness;
            
        end
        
        % Aerospike method
        %         function ctr = aerospike(thenozzle, variety)
        %               Coming Soon
        %         end
        
    end
    
end
