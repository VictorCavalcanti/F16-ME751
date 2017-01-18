function [engine] = simEngine3D_tsda(file)
%RUNLINE;
% [engine] = simEngine3D_tsda('tsda_test1');
%https://ocw.mit.edu/courses/mechanical-engineering/2-003j-dynamics-and-control-i-spring-2007/lecture-notes/lec21.pdf

% warning on verbose
warning off MATLAB:nearlySingularMatrix
acf = strcat(file,'.acf');
adm = strcat(file,'.adm');
%Run Dynamics analysis
tic
engine = simEngine3D(acf,adm);
engine = engine.runDynamics();
toc
plot(engine.myTimes,engine.myPosition(3,:),'b','LineWidth',1.5); hold on;
grid on;
xlabel('Time (s)');
ylabel('Displacement (m)');
title('simEngine3D Mass-Spring-Damper System');
legend('z-displacement of mass');
% m = 1;
% c = 1;
% k = 100;
% f = (-9.81/k);
% c1 = 0.5-f;
% % del = -9.81/k;
% % c1 = del-f;
% % wn = sqrt(4*m*k-c^2)/(2*m);
% wn = sqrt(k/m);
% % z = ((c1*exp((-c/(2*m)).*engine.myTimes).*cos(wn.*engine.myTimes)-f))-1.1962;
% z = (c1*exp((-c/(2*m)).*engine.myTimes).*cos(wn.*engine.myTimes))+f;
% 
% plot(engine.myTimes, z,'r'); 

m = 1;
c = 1;
k = 100;
f = (-9.81/k);
c1 = -0.5-f;
wn = sqrt(k/m);
z = (c1*exp((-c/(2*m)).*engine.myTimes).*cos(wn.*engine.myTimes))+f+1;
figure; hold on;
plot(engine.myTimes, z,'r','LineWidth',1.5); 
grid on;
xlabel('Time (s)');
ylabel('Displacement (m)');
title('Analytical Solution to Mass-Spring-Damper System');
legend('z-displacement of mass');
ts = length(engine.myITS);
for t = 1:ts
   error(t) = abs((z(t)-engine.myPosition(3,t))/z(t)).*100; 
end
figure; hold on;
plot(engine.myTimes, error(:), 'k');
ylabel('Error (%)');
xlabel('time (s)');
title('Error (%) of simEngine3D vs Analytical solution');
grid on;
%Plot x-y position of point B0
% ts = length(engine.myITS);
% plot(engine.myTimes,engine.myPosition(16,:),'b-'); hold on;
% plot(engine.myTimes,engine.myPosition(17,:),'r-');
% plot(engine.myTimes,engine.myPosition(15,:),'k-');
% P3pos = zeros(3,ts);
% for t = 1:ts
%     A = getA(engine.myPosition(11:14,t));
%     P3pos(:,t) = engine.myPosition(8:10,t)+A*[0.5;0;0];

% solution = xlsread('A02_solution.xlsx');
% errorx = (Bpos(1,:)'-solution(:,1))./solution(:,1);
% errory = (Bpos(2,:)'-solution(:,2))./solution(:,2);
% figure; hold on;
% plot(P3pos(1,:),'b-','Linewidth',1.5);
% plot(P3pos(2,:),'r-','Linewidth',1.5);
% figure; hold on;
% plot(errorx,'b--');
% plot(errory,'r--');
end