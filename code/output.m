function output(which, filename, x, y, curve)

if which == 1
    disp("Outputting "+length(x)+" sets of coordinates for ANSYS");
    
    fileID = fopen(filename+'.txt','w');
    for i = 1:length(x)
        fprintf(fileID, curve+" "+i+" "+x(i)+" "+y(i)+" 0\n");
    end
    fclose(fileID);
    
elseif which == 2
    disp("CSV")
elseif which == 3
    disp("Outputting "+length(x)+" sets of coordinates for SOLIDWORKS");
    
    fileID = fopen(filename+'.txt','w');
    for i = 1:length(x)
        fprintf(fileID,x(i)+" "+y(i)+" 0\n");
    end
    fclose(fileID);
else
    disp("INVALID FLAG")
end

end

