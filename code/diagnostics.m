% Rocket Nozzle Calculator - Diagnostics
%
% Written by Vinay Williams
% Started on 17/05/21

function diagnostics(nozzle)

fileID = fopen("properties.txt","w");
fprintf(fileID, "Force                 : "+ nozzle.force+"N\n");
fprintf(fileID, "Impulse               : "+ nozzle.impulse+"Ns\n");
fprintf(fileID, "Specific Impulse      : "+nozzle.specific_impulse+"s\n");
fprintf(fileID, "Exhaust Mass Flow Rate: "+nozzle.mass_flow_rate_exhaust+"kg/s\n");
fprintf(fileID, "Specific Heat Ratio   : "+nozzle.specific_heat_ratio+"\n");
fprintf(fileID, "Exhaust Molar Mass    : "+nozzle.molar_mass_exhaust+"kg\n");

fprintf(fileID, "Ambient Pressure      : "+nozzle.pressure_ambient+"Pa\n");

fprintf(fileID, "Chamber Pressure      : "+nozzle.pressure_chamber+"Pa\n");
fprintf(fileID, "Chamber Temperature   : "+nozzle.temperature_chamber+"K\n");

fprintf(fileID, "Burn Time             : "+nozzle.burn_time+"s\n");

fprintf(fileID, "Weigth Flow Rate      : "+nozzle.weight_flow_rate+"N/s\n");

fprintf(fileID, "Throat Radius         : "+nozzle.radius_throat+"m\n");
fprintf(fileID, "Throat Temperature    : "+nozzle.temperature_throat+"K\n");
fprintf(fileID, "Throat Area           : "+nozzle.area_throat+"m^2\n");
fprintf(fileID, "Throat Pressure       : "+nozzle.pressure_throat+"Pa\n");

fprintf(fileID, "Exit Area             : "+nozzle.area_exit+"m^2\n");
fprintf(fileID, "Exit Mach             : "+nozzle.exit_mach);
fprintf(fileID, "Exit Radius           : "+nozzle.radius_exit+"m\n");
fprintf(fileID, "Exit Speed of Sound   : "+nozzle.exit_speed_sound+"m/s\n");
fprintf(fileID, "Exit Mach             : "+nozzle.exit_mach+"\n");
fprintf(fileID, "Exit Velocity         : "+nozzle.velocity_exit+"m/s\n");
fprintf(fileID, "Exit Temperature      : "+nozzle.temperature_exit+"K\n");
fprintf(fileID, "Exit Pressure         : "+nozzle.pressure_exit+"Pa\n");

fclose(fileID);

disp("Force                 : "+ nozzle.force+"N")
disp("Impulse               : "+ nozzle.impulse+"Ns")
disp("Specific Impulse      : "+nozzle.specific_impulse+"s")
disp("Exhaust Mass Flow Rate: "+nozzle.mass_flow_rate_exhaust+"kg/s")
disp("Specific Heat Ratio   : "+nozzle.specific_heat_ratio)
disp("Exhaust Molar Mass    : "+nozzle.molar_mass_exhaust+"kg")

disp("Ambient Pressure      : "+nozzle.pressure_ambient+"Pa")

disp("Chamber Pressure      : "+nozzle.pressure_chamber+"Pa")
disp("Chamber Temperature   : "+nozzle.temperature_chamber+"K")

disp("Burn Time             : "+nozzle.burn_time+"s")

disp("Weigth Flow Rate      : "+nozzle.weight_flow_rate+"N/s")

disp("Throat Radius         : "+nozzle.radius_throat+"m")
disp("Throat Temperature    : "+nozzle.temperature_throat+"K")
disp("Throat Area           : "+nozzle.area_throat+"m^2")
disp("Throat Pressure       : "+nozzle.pressure_throat+"Pa")

disp("Exit Area             : "+nozzle.area_exit+"m^2")
disp("Exit Mach             : "+nozzle.exit_mach)
disp("Exit Radius           : "+nozzle.radius_exit+"m")
disp("Exit Speed of Sound   : "+nozzle.exit_speed_sound+"m/s")
disp("Exit Mach             : "+nozzle.exit_mach)
disp("Exit Velocity         : "+nozzle.velocity_exit+"m/s")
disp("Exit Temperature      : "+nozzle.temperature_exit+"K")
disp("Exit Pressure         : "+nozzle.pressure_exit+"Pa")


end

