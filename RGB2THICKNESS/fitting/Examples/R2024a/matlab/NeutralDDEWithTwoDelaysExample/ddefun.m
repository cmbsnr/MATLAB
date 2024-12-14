function yp = ddefun(t,y,ydel,ypdel) 
    yp = 1 + y - 2*ydel^2 - ypdel;
end