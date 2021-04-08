ccc

T = 4;
x = linspace(0,2*pi*T,100*T);


y =@(A,x) cos(A*x);

f = 2;

plot(x,y(f,x))
