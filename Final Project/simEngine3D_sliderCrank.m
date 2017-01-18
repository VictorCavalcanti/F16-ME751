function [engine] = simEngine3D_sliderCrank(file)
%RUNLINE;
% [engine] = simEngine3D_sliderCrank('sliderCrank');

%Bars have circular cross section, 1m length, 1kg weight, slender 
acf = strcat(file,'.acf');
adm = strcat(file,'.adm');
%Run Dynamics analysis
tic
engine = simEngine3D(acf,adm);
engine = engine.runDynamics();
toc

plot(engine.myTimes(:),engine.myPosition(15,:),'b-','LineWidth',1.5); hold on;
xlabel('Time (s)');
ylabel('Position (m)');
title('Slider-Crank');
legend('x-displacement of slider');
grid on;
%Plot x-y position of point B0
% ts = length(engine.myITS);
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