g = xlsread("inputs.xlsx","CD","C3");   % Gravitational Aceeleration
R = xlsread("inputs.xlsx","CD","C4");   % Universal Gas Constant

% Inputs
force = xlsread("inputs.xlsx","CD","C6");   % Force desired             [N]
impulse = xlsread("inputs.xlsx","CD","C7"); % Total Impulse             [Ns]
mdot = xlsread("inputs.xlsx","CD","C8");    % Mass Flow                 [kg/s]
gamma = xlsread("inputs.xlsx","CD","C9");   % Specific Heat Ratio
pc = xlsread("inputs.xlsx","CD","C10");     % Chamber Pressure          [bar]
tc = xlsread("inputs.xlsx","CD","C11");     % Chamber Temperature       [kelvin]
M = xlsread("inputs.xlsx","CD","C12");      % Molar Mass of Exhaust Gas [g]
alt = xlsread("inputs.xlsx","CD","C13");    % Operational Altitude      [km]