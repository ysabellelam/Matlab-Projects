function [T, X, Y, Z, U, V, W] = missile(X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn);

load('terrain.mat')

dt = .005;
m = m0;
r = .2;
A = pi*r^2;
g = 9.81;
d = 1.2; 

T = 0;
X = zeros(length(T)); X(1) = X0;
Y = zeros(length(T)); Y(1) = Y0;
Z = zeros(length(T)); Z(1) = Z0;
U = zeros(length(T)); U(1) = 0;
V = zeros(length(T)); V(1) = 0;
W = zeros(length(T)); W(1) = 0;

n = 1;

while Z(n) >= interp2(x_terrain,y_terrain,h_terrain,X(n),Y(n))
    [Th_x(n), Th_y(n), Th_z(n)] = thrust(T(n), Thmag0, theta, phi, Tburn, U(n), V(n), W(n));
    
    T(n+1) = T(n) +dt;
    m(n+1) = mass(T(n+1),m0,mf,Tburn);
    
    Vmag = sqrt(U(n)^2+V(n)^2+W(n)^2);
    Cd(n) = drag_coeff(Vmag);
    
    U(n+1) = U(n) + ((Th_x(n)/m(n)) - (((d*A)/(2*m(n)))*U(n)*Cd(n))*Vmag)*dt;
    V(n+1) = V(n) + ((Th_y(n)/m(n)) - (((d*A)/(2*m(n)))*V(n)*Cd(n))*Vmag)*dt;
    W(n+1) = W(n) + ((Th_z(n)/m(n)) - (((d*A)/(2*m(n)))*W(n)*Cd(n))*Vmag-g)*dt;
    
    X(n+1) = X(n) +U(n+1)*dt;
    Y(n+1) = Y(n) +V(n+1)*dt;
    Z(n+1) = Z(n) + W(n+1)*dt; 
    
    n = n+1;
    
end

end
    

