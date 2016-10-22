function [engine] = simEngine3D_A6P2(file)
%Assignment 6-3 - simEngine3D-A6P3
%%
%RUN LINE:     [engine] = simEngine3D_A6P2('revJoint');
%%
acf = strcat(file,'.acf');
mdl = strcat(file,'.mdl');
engine = simEngine3D(acf,mdl);
% engine.myFlags = [1, 1, 0, 0];
% engine.myPosition(1:14,1) = [engine.myBodies(1).q0'; 0;0;0;1;0;0;0];
% engine.myVelocity = zeros(14,1);
engine.myPosition(1:7,1) = engine.myBodies(1).q0;
engine = engine.getPhiJac(1);
disp('Phi');
disp(engine.myPhi);
disp('Jac');
disp(engine.myJac);

% engine = engine.PositionAnalysis();

end
