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
                        
plt = 1;                % 1   = Yes
                        % <1> = No
                        
plt_save = 1;           % 1   = Yes
                        % <1> = No
                        
contour_solver = 2;     % 1 = Method of Characteristics
                        % 2 = Rao (Not operational yet)
                        % 3 = Conical
                       
   
% Only for conical nozzle and rao
div_ang = 15;       % Diverging Angle           [deg]
conv_ang = 30;      % Converging Angle          [deg]
                        
% Only for Rao method
s = 0.8;

out_fmt = 3;            % 1 = ANSYS
                        % 2 = CSV
                        % 3 = SOLIDWORKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noz = nozzle(force, impulse, mdot, gamma, pc, M, tc, alt, tcknss, ...
                contour_solver, conv_ang, div_ang, s, "test");
            
points = noz.contour;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Console Diagnostics

noz.printdiagnostics(noz)
noz.writediagnostics(noz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output

noz.write(noz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
