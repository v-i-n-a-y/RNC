function pm = pm(me, g)
    A = sqrt((g+1)/(g-1));
    B = (g-1)/(g+1);
    pm =  A*atan(sqrt(B*(me^2-1))) - atan(sqrt(me^2-1));
end