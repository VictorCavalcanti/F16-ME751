function [engine] = simEngine3D_fourBar(file)
%RUNLINE;
% [engine] = simEngine3D_fourBar('fourBarN2');
% [engine] = simEngine3D_fourBar('fourBar');
%Bars have circular cross section, 1m length, 1kg weight, slender rods
%%%% NEED TO RELAX CONSTRAINTS. MAKING BOTH GROUND ATTACHMENTS SPHERICAL
%%%% WORKS OR PERHAPS 3 REV ONLY X & Y CD ON GROUND ON $ ALSO WORKS
acf = strcat(file,'.acf');
adm = strcat(file,'.adm');
%Run Dynamics analysis
tic
engine = simEngine3D(acf,adm);
engine = engine.runDynamics();
toc

%Plot x-y position of point B0
ts = length(engine.myITS);
Bpos = zeros(3,ts);
for t = 1:ts
    A = getA(engine.myPosition(4:7,t));
    Bpos(:,t) = engine.myPosition(1:3,t)+A*[0.5;0;0];
end
% solution = xlsread('A02_solution.xlsx');
% errorx = (Bpos(1,:)'-solution(:,1))./solution(:,1);
% errory = (Bpos(2,:)'-solution(:,2))./solution(:,2);
figure; hold on;
plot(engine.myTimes(1:end),Bpos(1,1:end),'r-','Linewidth',1.5);
plot(engine.myTimes(1:end),Bpos(2,1:end),'b--','Linewidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Displacement (m)');
title('Single 4-Bar Mechanism');
legend(['x-displacement of point B0'; 'y-displacement of point B0']);
end
% figure; hold on;
% plot(errorx,'b--');
% plot(errory,'r--');