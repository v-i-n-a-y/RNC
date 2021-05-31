% Rocket Nozzle Calculator - Main
%
% Written by Vinay Williams
% Started on 17/02/21
%
% See documentation for details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
clear; clc; home; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flags
xl = 1;                 % 1   = Get inputs from excel file (in work)
                        % <1> = Use variables defined in script
                        
plt = 1;                % 1   = Yes
                        % <1> = No
                        
plt_save = 1;                % 1   = Yes
                        % <1> = No
                        
contour_solver = 1;     % 1 = Method of Characteristics
                        % 2 = Rao (Not operational yet)
                        % 3 = Conical
                       
    if contour_solver == 3
        % Only for conical nozzle
        div_ang = 15;       % Diverging Angle           [deg]
        conv_ang = 30;      % Converging Angle          [deg]
    elseif contour_solver == 2 
        % Only for Rao method
        s = 0.8;
    end
    
out_fmt = 3;            % 1 = ANSYS
                        % 2 = CSV
                        % 3 = SOLIDWORKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
g = 9.81;               % Gravitational Constant [m/s/s]
R = 8.314;              % Universal Gas Constant [J/mol/k]

% Nozzle Performance Parameters
force = 8000;           % Force desired             [N]
impulse = 40000;        % Total Impulse             [Ns]
mdot = 2;               % Mass Flow                 [kg/s]
gamma = 1.2;            % Specific Heat Ratio
pc = 4;                 % Chamber Pressure          [Bar]
tc = 3200;              % Chamber Temperature       [K]
M = 27;                 % Molar Mass of Exhaust Gas [g]
alt = 10;              % Operational Altitude      [km]

% Nozzle Mechanical Parameters
tcknss = 2;             % Wall Thickness            [mm]

if xl == 0
    disp("NEEDS TO BE FINSHED/FIXED")
else
    disp("Using Parameters defined in script")
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overarching Calculations
calculations;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Contour 
if contour_solver == 1
    moc;
elseif contour_solver == 2
    rao;
elseif contour_solver == 3
    conical; 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Console Diagnostics
diagnostics(nozzle)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
output(out_fmt,"outer-points", nozzle.xpoints.outer,...
        nozzle.ypoints.outer, 1)
    
output(out_fmt,"inner-points", nozzle.xpoints.inner,...
        nozzle.ypoints.inner, 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
