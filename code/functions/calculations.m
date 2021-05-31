% Conversion
M = M/1000;                 % Conver to kilograms
pc = pc*100000;             % Convert to pascals
R = R/M;                    % Universal gas constant/molar mass
alt = alt *1000;            % Convert to meters
tcknss = tcknss/1000;       % Convert to millimeters

% Housekeeping
nozzle.force = force;
nozzle.impulse = impulse;
nozzle.mass_flow_rate_exhaust = mdot;
nozzle.specific_heat_ratio = gamma;
nozzle.pressure_chamber = pc;
nozzle.molar_mass_exhaust = M;
nozzle.temperature_chamber = tc;
nozzle.altitude = alt;

if contour_solver == 3
    nozzle.converging_angle = conv_ang;
    nozzle.diverging_angle = div_ang;
    clear conv_ang div_ang
elseif contour_solver == 2
    nozzle.ratio = s;
    clear s
end 

nozzle.wall_thickness = tcknss;
clear force impulse mdot tc M gamma pc alt conv_ang div_ang tcknss

% Calculations

if (11000>nozzle.altitude) && (nozzle.altitude<25000)
    nozzle.temperature_ambient = -56.46; %C
    nozzle.pressure_ambient = 1000*(22.65*exp(1.73-0.000157*nozzle.altitude));
elseif nozzle.altitude>=25000
    nozzle.temperature_ambient = -131.21 + 0.00299*nozzle.altitude ;
    nozzle.pressure_ambient = 1000*(2.488*((nozzle.temperature_ambient+273.1)/216.6)^-11.388);
else 
    nozzle.temperature_ambient = 15.04 - 0.00649*nozzle.altitude;
    nozzle.pressure_ambient = 1000*(101.29*((nozzle.temperature_ambient+273.1)/288.08)^5.256);
end

nozzle.burn_time = nozzle.impulse/nozzle.force; 

nozzle.specific_impulse = nozzle.force/nozzle.mass_flow_rate_exhaust; 

nozzle.weight_flow_rate = nozzle.mass_flow_rate_exhaust * g;

nozzle.temperature_throat   = tt(nozzle.temperature_chamber, ...
                                 nozzle.specific_heat_ratio);
                             
nozzle.pressure_throat      = pt(nozzle.pressure_chamber,...
                                 nozzle.specific_heat_ratio);
                             
nozzle.area_throat          = at(nozzle.mass_flow_rate_exhaust,...
                                 nozzle.pressure_throat, R, ...
                                 nozzle.temperature_throat, ...
                                 nozzle.specific_heat_ratio, g);
                             
nozzle.radius_throat        = sqrt(nozzle.area_throat/pi);

nozzle.exit_mach            = me(nozzle.specific_heat_ratio, ...
                                 nozzle.pressure_chamber, ...
                                 nozzle.pressure_ambient);
                             
nozzle.area_exit            = ae(nozzle.area_throat,nozzle.exit_mach,...
                                 nozzle.specific_heat_ratio);
                             
nozzle.radius_exit          = sqrt(nozzle.area_exit/pi);

nozzle.pressure_exit        = pe(nozzle.pressure_throat, ...
                                 nozzle.exit_mach, ...
                                 nozzle.specific_heat_ratio);
                             
nozzle.temperature_exit     = te(nozzle.temperature_throat, ...
                                 nozzle.exit_mach, ...
                                 nozzle.specific_heat_ratio);
                             
nozzle.exit_speed_sound     = a(nozzle.specific_heat_ratio,...
                                R ,nozzle.temperature_exit);
                            
nozzle.velocity_exit        = nozzle.exit_speed_sound * nozzle.exit_mach;

nozzle.pressure_ratio       = nozzle.pressure_ambient/...
                                nozzle.pressure_chamber;
                            
nozzle.temperature_ratio    = (2*nozzle.specific_heat_ratio*R*...
                                nozzle.temperature_chamber)/(...
                                nozzle.specific_heat_ratio-1);
                            
nozzle.pressure_ratio_2     = (nozzle.pressure_ambient/...
                                nozzle.pressure_chamber)^((...
                                nozzle.specific_heat_ratio-1)/...
                                nozzle.specific_heat_ratio);

nozzle.velocity_throat      = sqrt((2*nozzle.specific_heat_ratio*R*...
                                nozzle.temperature_chamber)/(...
                                nozzle.specific_heat_ratio+1));  
                            
nozzle.area_ratio           = nozzle.area_exit/nozzle.area_throat;