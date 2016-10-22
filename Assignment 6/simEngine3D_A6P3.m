function [engine] = simEngine3D_A6P3(file)
%Assignment 6-3 - simEngine3D-A6P3
%%
%RUN LINE: [engine] = simEngine3D_A6P3('revJoint')
%%
acf = strcat(file,'.acf');
mdl = strcat(file,'.mdl');
engine = simEngine3D(acf,mdl);
engine = engine.PositionAnalysis();
engine = engine.VelocityAnalysis();
engine = engine.AccelerationAnalysis();


%Gets information for point Q (the base of the pendulum)
figure; 
subplot(3,1,1);
plot(engine.myTimes(2:end),engine.myPosition(1,:),'b','LineWidth',1.5); 
hold on; 
plot(engine.myTimes(2:end),engine.myPosition(2,:),'r','LineWidth',1.5); 
plot(engine.myTimes(2:end),engine.myPosition(3,:),'k','LineWidth',1.5); 
legend('rx-1','ry-1','rz-1');
xlabel('time (seconds)');
ylabel('units');
title('Position - LRF- Oxyz_Prime') ;

subplot(3,1,2);
plot(engine.myTimes(2:end),engine.myVelocity(1,:),'b','LineWidth',1.5); 
hold on; 
plot(engine.myTimes(2:end),engine.myVelocity(2,:),'r','LineWidth',1.5);
plot(engine.myTimes(2:end),engine.myVelocity(3,:),'k','LineWidth',1.5);
legend('vx-1','vy-1','vz-1');
xlabel('time (s)');
ylabel('units/s');
title('Velocity - LRF- Oxyz_Prime');

subplot(3,1,3);
plot(engine.myTimes(2:end),engine.myAcceleration(1,:),'b','LineWidth',1.5);
hold on; 
plot(engine.myTimes(2:end),engine.myAcceleration(2,:),'r','LineWidth',1.5);
plot(engine.myTimes(2:end),engine.myAcceleration(3,:),'k','LineWidth',1.5);
legend('ax-1','ay-1','az-1');
xlabel('time (s)');
ylabel('units/s^2')
title('Acceleration - LRF- Oxyz_Prime') ;
%Get quantities for point Q in body 1
% sQi = ri+Ai*sQibar;
% sQidot = rid+B_p_sqi*pid;
% sQiddot = ridd+B_pd_sqi*pd+B_p_sqi*pdd;
sQiBAR = [-2,0,0]';
qPosition = zeros(3,1000);
qVelocity = zeros(3,1000);
qAcceleration = zeros(3,1000);
for i = 1:1000;
    r = engine.myPosition(1:3,i);
    p = engine.myPosition(4:7,i);
    rd = engine.myVelocity(1:3,i);
    pd = engine.myVelocity(4:7,i);
    rdd = engine.myAcceleration(1:3,i);
    pdd = engine.myAcceleration(4:7,i);
    ei = p(2:end);
    eiTIL = [ 0 -ei(3) ei(2); ei(3) 0 -ei(1); -ei(2) ei(1) 0];
    Ai = (p(1).^2-(ei')*ei)*eye(3)+2*ei*(ei')+2*p(1)*eiTIL;
    B_p_sqi = getB(p,sQiBAR);
    B_pd_sqi = getB(pd,sQiBAR);
    qPosition(:,i) = r+Ai*sQiBAR;
    qVelocity(:,i) = rd+B_p_sqi*pd;
    qAcceleration(:,i) = rdd+B_pd_sqi*pd+B_p_sqi*pdd;
end
% PLOT VALUES FOR POINT Q
figure; 
subplot(3,1,1);
plot(engine.myTimes(2:end),qPosition(1,:),'b','LineWidth',1.5); 
hold on; 
plot(engine.myTimes(2:end),qPosition(2,:),'r','LineWidth',1.5); 
plot(engine.myTimes(2:end),qPosition(3,:),'k','LineWidth',1.5); 
legend('rQx-1','rQy-1','rQz-1');
xlabel('time (s)');
ylabel('units');
title('Position - sQi') 

subplot(3,1,2);
plot(engine.myTimes(2:end),qVelocity(1,:),'b','LineWidth',1.5); 
hold on; 
plot(engine.myTimes(2:end),qVelocity(2,:),'r','LineWidth',1.5);
plot(engine.myTimes(2:end),qVelocity(3,:),'k','LineWidth',1.5);
legend('vQx-1','vQy-1','vQz-1');
xlabel('time (s)');
ylabel('units/s')
title('Velocity - sQi') 
subplot(3,1,3);
plot(engine.myTimes(2:end),qAcceleration(1,:),'b','LineWidth',1.5);
hold on; 
plot(engine.myTimes(2:end),qAcceleration(2,:),'r','LineWidth',1.5);
plot(engine.myTimes(2:end),qAcceleration(3,:),'k','LineWidth',1.5);
legend('aQx-1','aQy-1','aQz-1');
xlabel('time (s)');
ylabel('units/s^2')
title('Acceleration - sQi') 
end
