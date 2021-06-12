function ae = ae(at, me, y)
y1 = y-1;
y2 = y+1;
ae = (at/me) *((1+(y1/2)*me^2)/(y2/2))^(y2/(2*y1));
end

