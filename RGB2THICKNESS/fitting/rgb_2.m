wavelength = [600;560;470]/1.4;%RGB光波长
t = 1:500;
r = 20+80*exp(-2*0.0002*t)+80*cos(4*pi*t/wavelength(1)).*exp(-2*0.0002*t);
g = 20+80*exp(-2*0.0003*t)+80*cos(4*pi*t./wavelength(2)).*exp(-2*0.0002*t);
b = 20+80*exp(-2*0.0004*t)+80*cos(4*pi*t./wavelength(3)).*exp(-2*0.0002*t);
rgbmatrix = [r;g;b]