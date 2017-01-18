function [engine] = simEngine3D_simplePendulumn(file)
%RUNLINE;
% [engine] = simEngine3D_simplePendulumn('simplePendulumn');

acf = strcat(file,'.acf');
adm = strcat(file,'.adm');
%Run Dynamics analysis
tic
engine = simEngine3D(acf,adm);
engine = engine.runDynamics();
toc
%Solution from website is terrible looking, compare only first 10s with
%given graph of 10s of solution.
% soln = xlsread('A01_solution.xlsx');
% plot(soln(1:100,2),'b-'); %50s of simulation
plot(engine.myTimes(1:1000),engine.myPosition(1,1:1000),'r-','LineWidth',1.5);
hold on;
plot(engine.myTimes(1:1000),engine.myPosition(2,1:1000),'b--','LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Position (m)');
legend(['x-displacement of mass'; 'y-displacement of mass']);
title('Simple Pendulumn');
end