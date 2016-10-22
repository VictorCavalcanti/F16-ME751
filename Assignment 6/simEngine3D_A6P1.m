function [engine,engine2] = simEngine3D_A6P1(file)
%Assignment 6-1 - simEngine3D-A6P1
%%
%RUN LINE: [engine,engine2] = simEngine3D_A6P1('me751')
%%
% engine.myPhi ; engine.myJac; engine.myNu; engine.myGamma ~cons_dp2 
% engine2.myPhi ... ~ cons_d 


acf = strcat(file,'.acf');
mdl = strcat(file,'.mdl');
engine = simEngine3D(acf,mdl);
engine.myFlags = [1, 1, 1, 1];
engine.myPosition(1:14,1) = [engine.myBodies(1).q0'; 0;0;0;1;0;0;0];
engine.myVelocity = zeros(14,1);
engine2 = engine;
[engine.myPhi,engine.myJac,engine.myNu,engine.myGamma] = ...
                cons_dp2(engine.myKinCon{1},engine.myTimeStep,...
                engine.myFunTimes,engine.myPosition(:,1),...
                    engine.myVelocity(:,1),engine.myFlags);
%                 
%                 (obj.myKinCon{i},obj.myTimeStep,...
%                         obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags)
                
[engine2.myPhi,engine2.myJac,engine2.myNu,engine2.myGamma] = ...
cons_d(engine2.myKinCon{2},engine2.myTimeStep,engine.myFunTimes,...
engine2.myPosition(:,1),engine2.myVelocity(:,1),engine.myFlags);
end
