clear all; close all; clc; format long; 

name='Ysabelle Lam';
id='A13868700';
hw_num='project';

%% Task 1 Part 1
figure(1); hold on;
load terrain.mat;
surf(x_terrain/1000, y_terrain/1000, h_terrain/1000); shading interp; 

for m_id = 1:7
[X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn] = read_input('missile_data.txt', m_id);
[T, X, Y, Z, U, V, W] = missile(X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn);
load terrain.mat;
color = 'krbgmcy';
plot3(X/1000, Y/1000, Z/1000, color(m_id))
end
view(3);
legend('terrain', 'M 1','M 2','M 3','M 4','M 5','M 6','M 7');

for m_id = 1:7
[X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn] = read_input('missile_data.txt', m_id);
[T, X, Y, Z, U, V, W] = missile(X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn);
load terrain.mat;
color = 'krbgmcy';
plot3((X(end,end))/1000, (Y(end,end))/1000, (Z(end,end))/1000, [color(m_id) 'o'], 'MarkerSize', 5,'MarkerFaceColor', color(m_id), 'MarkerEdgeColor', color(m_id))
end

hold off;
xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)');
view(3); axis([0 30 0 30 0 3.5]); grid on; 
set(gca,'LineWidth',2,'FontSize',16);
title('Missile Trajectories');

%% Task 1 Part 2
dt = .005;
g= 9.81;
for m_id = 1:7
figure(2); hold on;
[X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn] = read_input('missile_data.txt', m_id);
[T, X, Y, Z, U, V, W] = missile(X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn);
color = 'krbgmcy';    
Vmag = sqrt(U.^2 + V.^2 + W.^2);
Ma = Vmag/340;

subplot(2,1,1)       
plot(T,Ma,color(m_id),'linewidth',2)
axis([0 60 0 8])
ylabel('Ma');
hold off;
end

title('Evolution of speed and acceleration');
grid on;

for m_id = 1:7
figure(2); hold on;
[X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn] = read_input('missile_data.txt', m_id);
[T, X, Y, Z, U, V, W] = missile(X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn);
color = 'krbgmcy';    
Vmag = sqrt(U.^2 + V.^2 + W.^2);

Acc = diff(Vmag)/dt;
Acc = Acc/g;

lAcc = length(Acc);

subplot(2,1,2)       
plot(T(:,1:lAcc), Acc,color(m_id),'linewidth',2)     
axis([0 60 -100 50])
xlabel('Time(s)'); ylabel('Acc/g');
hold off;
end
legend('M 1','M 2','M 3','M 4','M 5','M 6','M 7');
grid on; 

%% Task 1 Part 3
for m_id = 1:7
[X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn] = read_input('missile_data.txt', m_id);
[T, X, Y, Z, U, V, W] = missile(X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn);
Vmag = sqrt(U.^2 + V.^2 + W.^2);
Ma = Vmag/340;
color = 'krbgmcy';

figure(3); hold on;
plot (Ma, Z/1000,color(m_id))

end
hold off 
legend('M 1','M 2','M 3','M 4','M 5','M 6','M 7');

for m_id = 1:7
[X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn] = read_input('missile_data.txt', m_id);
[T, X, Y, Z, U, V, W] = missile(X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn);
Vmag = sqrt(U.^2 + V.^2 + W.^2);
Ma = Vmag/340;
color = 'krbgmcy';
figure(3); hold on; 

F= find(Ma>=1 & Ma<= 1.001);
H = F(1,1);
G = F(1,end);

plot(1, Z(H)/1000, [color(m_id) 'o'], 'MarkerSize', 5,'MarkerFaceColor', color(m_id), 'MarkerEdgeColor', color(m_id))
plot(1, Z(G)/1000, [color(m_id) 'o'], 'MarkerSize', 5,'MarkerFaceColor', color(m_id), 'MarkerEdgeColor', color(m_id))

end
hold off;
axis([0 9 0 3.5])
grid on;
title('Sonic Barrier');

%% Task 2

dt = .005;
g= 9.81;
for m_id = 1:7
[X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn] = read_input('missile_data.txt', m_id);
[T, X, Y, Z, U, V, W] = missile(X0, Y0, Z0, m0, mf, Thmag0, theta, phi, Tburn);
Vmag = sqrt(U.^2 + V.^2 + W.^2);
Ma = Vmag/340;
Acc = diff(Vmag)/dt;
Acc = Acc/g;
    missile_ID(m_id) = {m_id};
    landing_time(m_id)={T(1,end)};
    travel_distance(m_id)= {sum(sqrt(diff(U).^2 + diff(V).^2 + diff(W).^2 ))};
    max_height_position(m_id) = {[max(X), max(Y), max(Z)]};
    max_height_Ma(m_id)={Ma(find(max(Z)==Z))};
    max_height_Acc(m_id) = {Acc(find(max(Z)==Z))};
    landing_location(m_id) = {[X(1,end), Y(1,end), Z(1,end)]}; 
    landing_Ma(m_id)= {Ma(1,end)};
    landing_Acc(m_id)= {Acc(1,end)};
    vmag(m_id) = {Vmag(1,end)};
    
end

flight_stat = struct('missile_ID',missile_ID,'landing_time',landing_time,...
    'travel_distance',travel_distance, 'max_height_position', max_height_position,...
    'max_height_Ma',max_height_Ma, 'max_height_Acc', max_height_Acc, ...
    'landing_location',landing_location, 'landing_Ma', landing_Ma,...
    'landing_Acc', landing_Acc);

%% Task 3
fileID = fopen('report.txt','w');

fprintf(fileID,'%s\n',name);
fprintf(fileID,'%s\n',id);
fprintf(fileID,'%s\n','M ID, landing time (s), travel distance (m), landing speed (m/s),landing acceleration (m/sˆ2)');
fprintf(fileID, '%d\t %15.9e\t %15.9e\t %15.9e\t %15.9e\n',missile_ID{1} , landing_time{1},travel_distance{1},vmag{1},landing_Acc{1}); 
fprintf(fileID, '%d\t %15.9e\t %15.9e\t %15.9e\t %15.9e\n',missile_ID{2} , landing_time{2},travel_distance{2},vmag{2},landing_Acc{2});
fprintf(fileID, '%d\t %15.9e\t %15.9e\t %15.9e\t %15.9e\n',missile_ID{3} , landing_time{3},travel_distance{3},vmag{3},landing_Acc{3});
fprintf(fileID, '%d\t %15.9e\t %15.9e\t %15.9e\t %15.9e\n',missile_ID{4} , landing_time{4},travel_distance{4},vmag{4},landing_Acc{4});
fprintf(fileID, '%d\t %15.9e\t %15.9e\t %15.9e\t %15.9e\n',missile_ID{5} , landing_time{5},travel_distance{5},vmag{5},landing_Acc{5});
fprintf(fileID, '%d\t %15.9e\t %15.9e\t %15.9e\t %15.9e\n',missile_ID{6} , landing_time{6},travel_distance{6},vmag{6},landing_Acc{6});
fprintf(fileID, '%d\t %15.9e\t %15.9e\t %15.9e\t %15.9e\n',missile_ID{7} , landing_time{7},travel_distance{7},vmag{7},landing_Acc{7});
fclose(fileID);
type report.txt;

%% END
p1a ='See figure 1';
p1b ='See figure 2';
p1c ='See figure 3';
p2a = flight_stat(1);
p2b = flight_stat(2);
p2c = flight_stat(3);
p2d = flight_stat(4);
p2e = flight_stat(5);
p2f = flight_stat(6);
p2g = flight_stat(7);
p3 = evalc('type report.txt');
