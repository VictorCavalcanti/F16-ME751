function [engine] = simEngine3D_A6P3(file)
%Assignment 6-3 - simEngine3D-A6P3
%[engine] = simEngine3D_A6P2('revJoint')
acf = strcat(file,'.acf');
mdl = strcat(file,'.mdl');
engine = simEngine3D(acf,mdl);
engine = engine.PositionAnalysis();
engine = engine.VelocityAnalysis();
engine = engine.AccelerationAnalysis();

%Plotting goes here, but I've been stuck with Jacobian of PHI being
%singular for many hours and I dont know how to get past it :(



end
