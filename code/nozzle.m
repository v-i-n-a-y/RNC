% RNC - ROCKET NOZZLE CALCULATOR
%
% Nozzle Class

classdef nozzle
    
    % Initally defined variables
    properties (Dependent = false)
        force
        impulse
        propellant_flow_rate
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
        name
        width                   % Linear Aerospike width
        
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
        diagnostics
        exit_density
        throat_density
        throat_mach
        throat_speed_sound
        exit_mass_flow
        throat_mass_flow
        characteristic_exhaust_velocity
    end
    
    properties (Constant, Access = private)
        g = 9.812;                % Gravitational Acceleration
        R = 8.314;                % Universal Gas Constant
        bar = 100000;             % Bar to pascal
        rtod = 180/pi;            % Radians to degrees
        dtor = pi/180;            % Degrees to radians
    end
   
    
    methods (Static)
        
        % Prandtl Meyer Function
        
        
        function out = SetGet_static(arg)
            
            persistent var
            
            if isempty(var)
                var = 0; 
            end
            if nargin == 1
                var = var+1;
            else
                
                out = var;
            end 
            
        end
    end
    
    methods
        
        % Method of characteristics method
        function ctr = moc(obj)
            
            % Prandtl Meyer Function
            A = sqrt((obj.y+1)/(obj.y-1));
            B = (obj.y-1)/(obj.y+1);
            v_PM = @(x) A*atan(sqrt(B*(x^2-1))) - atan(sqrt(x^2-1));
            
            % Diverging Angle
            theta_max = 0.5*v_PM(obj.exit_mach)*obj.rtod; % Max Angle
            delta_theta = (90-theta_max) - fix(90-theta_max); %
            n = theta_max*2;
            
            [ambienttemp, M, P, RR, LR, SL, x, y0, s, b, xw, yw] = deal(0);
            
            for m = 2:n+1
                
                ambienttemp(m) = (delta_theta + (m-1))*obj.dtor;
                
                x_int = [1 1.0001*obj.exit_mach];
                func = @(x) ambienttemp(m) - v_PM(x);
                M(m) = fzero(func,x_int);
                P(m) = 0 + obj.throat_radius*tan(ambienttemp(m)); %X-AXIS POINT
                RR(m) = -obj.throat_radius/P(m);
                LR(m) = tan(ambienttemp(m)+asin(1/M(m)));
                SL(m) = -RR(m);
            end
            
            P(1) = [];
            l = length(P);
            for j = 1:l
                P1 = [0 obj.throat_radius];
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
                x(c) = (obj.throat_radius+SL(c)*P(c))/(SL(c)-F);
                y0(c) = F*x(c)+obj.throat_radius;
                X_P = [P(c) x(c)];
                Y_P = [0 y0(c)];
                
%                 plot(X_P,Y_P,'b');
%                 plot(X_P,-Y_P,'b');
%                 hold on
               
            end
            
            ambient_temperatureM = theta_max*obj.dtor;
            
            xw(1) = (obj.throat_radius+SL(1)*P(1))/(SL(1)-tan(ambient_temperatureM));
            yw(1) = tan(ambient_temperatureM)*xw(1)+obj.throat_radius;
            X_P2 = [P(1) xw];
            Y_P2 = [P(2) yw];
            
%             plot(X_P2,Y_P2,'g');
%             plot(X_P2,-Y_P2,'g');
            
            
            %DIVIDE (delta slopes)
            delta_thetaW = tan(ambient_temperatureM)/(length(P)-1);
            s(1) = tan(ambient_temperatureM);
            b(1) = obj.throat_radius;
            
            
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
            yw = [obj.throat_radius yw];
            
            
%             plot(xw, yw, 'cyan')
%             plot(xw, -yw, 'cyan')
%             h = plot(NaN,NaN,'ocyan');
%             legend(h,'wall');
%             title("Bell Nozzle - Method of Characteristics")
           
            ctr.x.outer = xw;
            ctr.y.outer = yw+obj.wall_thickness;
            
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
        
        % Static function to plot the nozzle contour
        function plot(obj)
            close all
            hold on
            
            plot(obj.contour.x.inner, obj.contour.y.inner);
            plot(obj.contour.x.inner, -obj.contour.y.inner);
        end
        
        % Static function to write contour points
        function write(obj, format, curve)
            contour1 = obj.contour; 
           
            if format == 1
                disp("Outputting "+length(contour1.x.inner)+" sets of coordinates for ANSYS");
                
                fileID = fopen(obj.name+'.txt','w');
                for i = 1:length(contour1.x.inner)
                    fprintf(fileID, curve+" "+i+" "+contour1.x.inner(i)+" "+contour1.y.inner(i)+" 0\n");
                end
                fclose(fileID);
                
            elseif format == 2
                disp("CSV")
            elseif format == 3
                disp("Outputting "+length(contour1.x.inner)+" sets of coordinates for SOLIDWORKS");
                
                fileID = fopen(obj.name+'.txt','w');
                for i = 1:length(contour1.x.inner)
                    fprintf(fileID,contour1.x.inner(i)+" "+contour1.y.inner(i)+" 0\n");
                end
                fclose(fileID);
            else
                disp("INVALID FLAG")
            end
        end
        
        % Function to write the diagnotics information on the nozzle to
        % file
        function writediagnostics(thenozzle)
            fileID = fopen(thenozzle.name+"properties.txt","w");
            fprintf(fileID, thenozzle.diagnostics);
            fclose(fileID);
        end
        
        % Static function to print diagnositcs to console
        function printdiagnostics(thenozzle)
            fprintf(thenozzle.diagnostics);
        end 
        
    end
   
    methods
        
        % Constructor
        function thisnozzle = nozzle(force, impulse, ...
                propellant_flow_rate, y, chamber_pressure, ...
                molar_mass_exhaust,chamber_temperature,altitude, ...
                wall_thickness,kind, converging_angle, diverging_angle, ...
                ratio, name)
            
            thisnozzle.force = force;
            thisnozzle.impulse = impulse;
            thisnozzle.chamber_temperature = chamber_temperature;
            thisnozzle.altitude = altitude * 1000;
            thisnozzle.kind = kind;
            thisnozzle.molar_mass_exhaust = molar_mass_exhaust/1000;
            thisnozzle.propellant_flow_rate = propellant_flow_rate;
            thisnozzle.wall_thickness = wall_thickness;
            thisnozzle.chamber_pressure = chamber_pressure*100000;
            thisnozzle.y = y;
            thisnozzle.name = name;
            
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
            elseif kind == 4
%                 if ~exist('width','var')
%                     error("Linear Aerospike chosen: Please provide width")
%                 end
                thisnozzle.width = 10;
            end
            
            nozzle.SetGet_static(1)
            
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
            specific_impulse = thenozzle.force/thenozzle.propellant_flow_rate;
        end
        
        % Weight Flow Rate
        function weight_flow_rate = get.weight_flow_rate(thenozzle)
            weight_flow_rate = thenozzle.propellant_flow_rate * thenozzle.g;
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
            throat_area = (thenozzle.propellant_flow_rate/thenozzle.throat_pressure)*sqrt((thenozzle.R_dot*thenozzle.throat_temperature)/(thenozzle.y*thenozzle.g));
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
        
        % Throat Speed of Sound
        function throat_speed_sound = get.throat_speed_sound(thenozzle)
            throat_speed_sound = sqrt(thenozzle.y*thenozzle.R_dot*thenozzle.throat_temperature);
        end
        
        % Throat Mach
        function throat_mach = get.throat_mach(thenozzle)
            throat_mach = thenozzle.throat_velocity / thenozzle.throat_speed_sound;
        end
        
        % Area Ratio
        function area_ratio = get.area_ratio(thenozzle)
            %area_ratio = thenozzle.exit_area/thenozzle.throat_area;
            area_ratio = 1/thenozzle.exit_mach * ( ( (2/(thenozzle.y+1))...
                         * (1 + (thenozzle.y-1)/2 * thenozzle.exit_mach^2))...
                         ^ ((thenozzle.y+1)/(2*(thenozzle.y-1))) ) ;
        end
        
        % Characteristic Exhaust Velocity
        function characteristic_exhaust_velocity = get.characteristic_exhaust_velocity(thenozzle)
            
            characteristic_exhaust_velocity = thenozzle.chamber_pressure*thenozzle.throat_area/thenozzle.exit_mass_flow;
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
                contour = aerospike(thenozzle, 1);
            elseif thenozzle.kind == 5
                disp("Aerospike: Linear")
                contour = aerospike(thenozzle, 2);
            else
                disp("Wrong Input")
            end
            
        end
       
        % Throat Density
        function throat_density = get.throat_density(thenozzle)
            throat_density = thenozzle.throat_pressure/(thenozzle.R_dot*thenozzle.throat_temperature);
        end
        
        % Exit Density
        function exit_density = get.exit_density(thenozzle)
            exit_density = thenozzle.exit_pressure/(thenozzle.R_dot*thenozzle.exit_temperature);
        end
        
        % Mass Flow rate at throat
        function throat_mass_flow = get.throat_mass_flow(thenozzle)
            throat_mass_flow = thenozzle.throat_density*thenozzle.throat_area*thenozzle.throat_velocity;
        end
        
        % Mass Flow Rate at exit
        function exit_mass_flow = get.exit_mass_flow(thenozzle)
            exit_mass_flow = thenozzle.exit_density*thenozzle.exit_area*thenozzle.exit_velocity;
        end
        
%         Aerospike method
        function ctr = aerospike(thenozzle, variety)
            
            A = sqrt((thenozzle.y+1)/(thenozzle.y-1));
            B = (thenozzle.y-1)/(thenozzle.y+1);
            v_PM = @(x) A*atan(sqrt(B*(x^2-1))) - atan(sqrt(x^2-1));
            
            N = 100;
            % Zero fill matrices
            M_x     = zeros(1,N);
            nu_x    = zeros(1,N);
            mu_x    = zeros(1,N);
            phi_x   = zeros(1,N);
            RxRe    = zeros(1,N);
            XxRe    = zeros(1,N);
            disp("test")
            for i = 1:N
                
                %Increment Mach number
                M_x(i)   = thenozzle.throat_mach + (i-1) * (thenozzle.exit_mach-1)/(N-1);     % Mach number
                
                % Calculate Mach wave parameters
                nu_x(i)     = v_PM(M_x(i));                   % Prandtl-Meyer angle [rad]
                mu_x(i)     = asin(1/M_x(i));              % Mach angle [rad]
                phi_x(i)    = v_PM(thenozzle.exit_mach) - nu_x(i) + mu_x(i);              % Angle of wave relative to axis [rad]
                
                % Calculate radial position of contour point
                a = (1 - 1/thenozzle.area_ratio * ( ...
                    ( (2/(thenozzle.y+1)) * (1+ (thenozzle.y - 1)/2 * M_x(i)^2) )^( (thenozzle.y+1)/(2*(thenozzle.y-1)) ) ...
                    * sin(phi_x(i))) );
                
                if variety == 1 % Radial case
                    
                    RxRe(i) = sqrt(a);     % Non-dimensional radial position
                    XxRe(i)  = (1 - RxRe(i)) / tan(phi_x(i));   % Non-dimensional axial position
                    
                elseif variety == 2 % Linear case
                    
                    RxRe(i) = a;           % Non-dimensional radial position
                    XxRe(i)  = (1 - RxRe(i)) / tan(phi_x(i));   % Non-dimensional axial position
                    
                elseif variety == 3 % Linear case
                    
                    RxRe(i) = 1 - i/N;
                    XxRe(i) = 2.779 * i/N;
                end
                
                
            end
                         
                
                if variety == 1     % Radial case
                    
                    A_e = pi*(thenozzle.exit_radius)^2;  
                    
                elseif variety == 2 % Linear case
                    
                    A_e = 2*thenozzle.width*thenozzle.exit_radius;     
                    
                end
                
                % THROAT CONDITIONS
                P_t     = thenozzle.throat_pressure;
                v_t     = thenozzle.throat_velocity;
                delta   = pi/2 - v_PM(thenozzle.exit_mach);                          % Angle of flow relative to lip plane [rad]
                
                % CONDITIONS ALONG CONTOUR
                
                % Zero fill matrices
                M_x     = zeros(1,N);
                nu_x    = zeros(1,N);
                mu_x    = zeros(1,N);
                phi_x   = zeros(1,N);
                theta_x = zeros(1,N);
                R_x     = zeros(1,N);
                X_x     = zeros(1,N);
                A_x     = zeros(1,N);
                P_x     = zeros(1,N);
                T_x     = zeros(1,N);
                v_x     = zeros(1,N);
                ux_x    = zeros(1,N);
                uy_x    = zeros(1,N);
                F_x     = zeros(1,N);
                Isp     = zeros(1,N);
                
                % Contour loop
                for i = 1:N
                    
                    %Increment Mach number
                    M_x(i)      = thenozzle.throat_mach + (i-1) * (thenozzle.exit_mach-1)/(N-1);          % Mach number
                    
                    % Calculate Mach wave parameters
                    nu_x(i)     = v_PM(M_x(i));                        % Prandtl-Meyer angle [rad]
                    mu_x(i)     = asin(1/M_x(i));                   % Mach angle [rad]
                    phi_x(i)    = v_PM(thenozzle.exit_mach) - nu_x(i) + mu_x(i);             % Angle of wave relative to axis [rad]
                    theta_x(i)  = delta + nu_x(i);                      % Flow angle relative to nozzle axis [rad]
                    
                    %Calculate flow conditions
                    P_x(i)  = thenozzle.chamber_pressure/(1 + (thenozzle.y-1)/2 * M_x(i)^2)^(thenozzle.y/(thenozzle.y-1));   % Pressure [Pa]
                    T_x(i)  = thenozzle.chamber_temperature/(1 + (thenozzle.y-1)/2 * M_x(i)^2);             % Temperature [K]
                    v_x(i)  = M_x(i) * sqrt(thenozzle.y*thenozzle.R_dot*T_x(i));                % Velocity [m/s]
                    ux_x(i) = sin(theta_x(i))*v_x(i);                   % Velocity component in X direction [m/s]
                    uy_x(i) = -cos(theta_x(i))*v_x(i);                  % Velocity component in Y direction [m/s]
                    
                    %Calculate coordinates
                    R_x(i) = RxRe(i)*thenozzle.exit_radius;                               % Radial position [cm]
                    X_x(i) = XxRe(i)*thenozzle.exit_radius;                               % Axial position [cm]
                    
                    if variety == 1     % Radial case
                        
                        A_x(i)  = pi*R_x(i)^2;          % Spike cross sectional area [cm^2]
                        
                    elseif variety == 2 % Linear case
                        
                        A_x(i)  = 2*thenozzle.width*R_x(i);           % Spike cross sectional area [cm^2]
                        
                    end
                    
                    if i > 1
                        
                        if i == 2
                            
                            if variety == 1
                                
                                A_t = pi*(thenozzle.exit_radius^2-R_x(1)^2)/sin(phi_x(1));       % Throat area [mm^2]
                                
                            elseif variety == 2
                                
                                A_t = 2*thenozzle.width*(thenozzle.exit_radius-R_x(1))/sin(phi_x(1));          % Throat area [mm^2]
                                
                            end
                            
                            area_ratio2 = A_e/A_t;   % Actual area ratio
                            
                        end
                        
                        % Thrust contribution [N]
                        F_x(i) = (P_x(i-1) + P_x(i) - 2*thenozzle.ambient_pressure)/2 *(A_x(i-1)-A_x(i))/1e6;
                        
                        % Specific impulse [s]
                        Isp(i) = Isp(i-1) + ( ...
                            v_t/thenozzle.y ...
                            * 1/2 ...
                            * ((thenozzle.y+1)/2)^(thenozzle.y/(thenozzle.y-1)) ...
                            * ((P_x(i-1) + P_x(i) - 2*thenozzle.ambient_pressure)/thenozzle.chamber_pressure) ...
                            * (A_x(i-1)-A_x(i))/A_t...
                            * 1/thenozzle.g ...
                            );
                    end
                    
                end
                
                % ENGINE PERFORMANCE
                
                % Mass flow rate
                % If throat and exit flow rates are different, flow is likely not choked at throat
                m_dot_e = thenozzle.exit_mass_flow;
                m_dot   = m_dot_e;                              % Mass flow rate [kg/s]
                
                F       = m_dot*thenozzle.exit_velocity + (thenozzle.exit_pressure-thenozzle.ambient_pressure)*(A_e/1e6);      % Total thrust [N]
                
                F2      = sum(F_x) +  m_dot * v_t * sin(delta) ... % Total thrust [N]
                    + (P_t-thenozzle.ambient_pressure) * A_t/1e6 * sin(delta);
                
                Isp     = Isp + 1/thenozzle.g * v_t * sin(delta) ...            % Specific impusle [s]
                    * (1 + 1/thenozzle.y*(1 - ((thenozzle.y+1)/2)^(thenozzle.y/(thenozzle.y-1))*(thenozzle.ambient_pressure/thenozzle.chamber_pressure)) );
                
                Isp_F   = F/(m_dot*thenozzle.g);
                
                h_t = thenozzle.exit_radius *(thenozzle.area_ratio-sqrt(thenozzle.area_ratio-sin(delta)))/(thenozzle.area_ratio*sin(delta)); % Width of Throat
                
                
                % TRUNCATION
                
                % Find index of point where Isp at 90%
                i_tr    = find(Isp/Isp(N) > 0.90, 1);           % Index of truncation point
                
%                 % Find coordinate of truncation point
%                 X_tr    = X_x(i_tr);                            % Axial position of truncation [mm]
%                 XtrRe   = X_tr/thenozzle.exit_radius;                             % Non-dimensional axial position
%                 R_tr    = R_x(i_tr);                            % Radial position of truncation [mm]
%                 RtrRe   = R_tr/thenozzle.exit_radius;                             % Non-dimensional radial position
%                 
%                 % Find flow conditons at truncation point
%                 thenozzle.throat_machr    = M_x(i_tr);                            % Mach number at truncation point
%                 P_tr    = P_x(i_tr);                            % Pressure [s] at truncation point
%                 T_tr    = T_x(i_tr);                            % Temperature [s] at truncation point
%                 % Find truncated spike performance
%                 Isp_tr  = Isp(i_tr);                            % Specific impulse [s] of truncated spike
%                 F_tr    = Isp_tr*m_dot*g;                       % Thrust [N] of truncated spike
%                 
%                 
%                 % Pthenozzle.area_ratioT MASS ESTIMATION
%                 % Full spike
%                 L_fl    = max(X_x) - min(X_x);                  % Spike length [mm^3]
%                 V_fl    = trapz(X_x,A_x);                       % Spike volume [mm^3]
%                 m_fl    = V_fl*dens_s/(1e3);                    % Spike mass [g]
%                 
%                 % Truncated spike
%                 L_tr    = X_tr - min(X_x);                      % Spike length [mm^3]
%                 V_tr    = trapz(X_x(1:i_tr), A_x(1:i_tr));      % Spike volume [mm^3]
%                 m_tr    = V_tr*dens_s/(1e3);                    % Spike mass [g]
                
                ctr.x.inner = X_x;
                ctr.y.inner = R_x;
                
                

        end

        
        
        function diagnostics = get.diagnostics(thenozzle)
            diagnostics = ...
            "Force                 : "+thenozzle.force+"N\n"+...
            "Impulse               : "+thenozzle.impulse+"Ns\n"+...
            "Specific Impulse      : "+thenozzle.specific_impulse+"s\n"+...
            "Exhaust Mass Flow Rate: "+thenozzle.propellant_flow_rate+"kg/s\n"+...
            "Specific Heat Ratio   : "+thenozzle.y+"\n"+...
            "Exhaust Molar Mass    : "+thenozzle.molar_mass_exhaust+"kg\n"+...
            "Ambient Pressure      : "+thenozzle.ambient_pressure+"Pa\n"+...
            "Chamber Pressure      : "+thenozzle.chamber_pressure+"Pa\n"+...
            "Chamber Temperature   : "+thenozzle.chamber_temperature+"K\n"+...
            "Burn Time             : "+thenozzle.burn_time+"s\n"+...
            "Weigth Flow Rate      : "+thenozzle.weight_flow_rate+"N/s\n"+...
            "Throat Radius         : "+thenozzle.throat_radius+"m\n"+...
            "Throat Temperature    : "+thenozzle.throat_temperature+"K\n"+...
            "Throat Area           : "+thenozzle.throat_area +"m^2\n"+...
            "Throat Pressure       : "+thenozzle.throat_pressure+"Pa\n"+...
            "Exit Area             : "+thenozzle.exit_area+"m^2\n"+...
            "Exit Mach             : "+thenozzle.exit_mach+"\n"+...
            "Exit Radius           : "+thenozzle.exit_radius+"m\n"+...
            "Exit Speed of Sound   : "+thenozzle.exit_speed_sound+"m/s\n"+...
            "Exit Mach             : "+thenozzle.exit_mach+"\n"+...
            "Exit Velocity         : "+thenozzle.exit_velocity+"m/s\n"+...
            "Exit Temperature      : "+thenozzle.exit_temperature+"K\n"+...
            "Exit Pressure         : "+thenozzle.exit_pressure+"Pa\n";
        end
    end
    
end
