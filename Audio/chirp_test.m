t = 0:1/44100:5;
y = chirp(t,0,5,5000);
sound(y,44100)