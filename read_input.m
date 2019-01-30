function[X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn] = read_input(input_filename, M_id)
%This function reads the parameters stored in a given file
%inputs filename and missile id
%outputs initial position, mass, and thrust magnitude of missile as well as the final mass, angles theta and phi, and the fuel burn time 
missileinfo = importdata(input_filename, '\t', 7);

[r,c] = size([missileinfo.data]);
 

if M_id <= 0 || M_id >r
    X0 = NaN;
    Y0 =NaN;
    Z0 = NaN;
    m0 =NaN;
    mf= NaN;
    Thmag0 = NaN;
    theta = NaN;
    phi = NaN;
    Tburn = NaN;
    disp('Error: invalid Missile ID');
    
else 
    X0 = missileinfo.data(M_id,2);
    Y0 = missileinfo.data(M_id,3);
    Z0= missileinfo.data(M_id,4);
    m0= missileinfo.data(M_id,5);
    mf= missileinfo.data(M_id,6);
    Thmag0= missileinfo.data(M_id,7);
    theta= missileinfo.data(M_id,8);
    phi= missileinfo.data(M_id,9);
    Tburn= missileinfo.data(M_id,10);
end
