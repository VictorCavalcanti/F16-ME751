function [engine] = simEngine3D_A8P2(file)
% function [engine,torque] = simEngine3D_A8P2(file)
%RUNLINE;
% [engine] = simEngine3D_A8P2('doublePendulum');


acf = strcat(file,'.acf');
mdl = strcat(file,'.mdl');
engine = simEngine3D(acf,mdl);
% torque = zeros(3,numel(engine.myTimes-1));
engine = engine.runDynamics();
% torque(:,i) = 0.5*G*engine.myFreact(4:7,i);

% for i = 1:numel(engine.myTimes)-1
% p = engine.myPosition(4:7,i);
% e = p(2:end);
% eTIL = [ 0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];
% G = [-e, -eTIL+p(1)*eye(3)];
% A = (p(1).^2-e'*e)*eye(3)+2*e*e'+2*p(1).*eTIL;
% spi = engine.myKinCon{1}.sPiBAR;
% sPiTIL = [ 0 -spi(3) spi(2); spi(3) 0 -spi(1); -spi(2) spi(1) 0];
% % torque(:,i) = -0.5*G*engine.myFreact(4:7,i)-sPiTIL*A*engine.myFreact(1:3,i);
% torque(:,i) = 0.5*G*engine.myFreact(7:10,i);
% % torque(:,i) = A*torque(:,i);
% end
% plot(engine.myTimes(1:end),torque(1,:),'r','LineWidth',1.5); hold on;
% plot(engine.myTimes(1:end),torque(2,:),'k','LineWidth',1.5); 
% plot(engine.myTimes(1:end),torque(3,:),'b','LineWidth',1.5); 
% legend('Tx-1','Ty-1','Tz-1');
% xlabel('time (seconds)');
% ylabel('Torque (N*m)');
% title('Torque of Pendulum vs Time') ;


figure;
subplot(3,1,1); hold on;
plot(engine.myTimes(1:end),engine.myPosition(1,:),'r','LineWidth',1.5); 
plot(engine.myTimes(1:end),engine.myPosition(2,:),'k','LineWidth',1.5); 
plot(engine.myTimes(1:end),engine.myPosition(3,:),'b','LineWidth',1.5);
legend('Posx-1','Posy-1','Posz-1');
xlabel('time (seconds)');
ylabel('Position (m)');
title('Position of Body1 vs Time') ;
subplot(3,1,2); hold on;
plot(engine.myTimes(1:end),engine.myVelocity(1,:),'r','LineWidth',1.5); 
plot(engine.myTimes(1:end),engine.myVelocity(2,:),'k','LineWidth',1.5); 
plot(engine.myTimes(1:end),engine.myVelocity(3,:),'b','LineWidth',1.5);
legend('Velx-1','Vely-1','Velz-1');
xlabel('time (seconds)');
ylabel('Velocity (m/s)');
title('Velocity of Body1 vs Time') ;
subplot(3,1,3); hold on;
plot(engine.myTimes(1:end),engine.myAcceleration(1,:),'r','LineWidth',1.5); 
plot(engine.myTimes(1:end),engine.myAcceleration(2,:),'k','LineWidth',1.5); 
plot(engine.myTimes(1:end),engine.myAcceleration(3,:),'b','LineWidth',1.5);
legend('Accelx-1','Accely-1','Accelz-1');
xlabel('time (seconds)');
ylabel('Acceleration (m/s^2)');
title('Acceleration of Body1 vs Time') ;

figure;
subplot(3,1,1); hold on;
plot(engine.myTimes(1:end),engine.myPosition(8,:),'r','LineWidth',1.5); 
plot(engine.myTimes(1:end),engine.myPosition(9,:),'k','LineWidth',1.5); 
plot(engine.myTimes(1:end),engine.myPosition(10,:),'b','LineWidth',1.5);
legend('Posx-2','Posy-2','Posz-2');
xlabel('time (seconds)');
ylabel('Position (m)');
title('Position of Body2 vs Time') ;
subplot(3,1,2);  hold on;
plot(engine.myTimes(1:end),engine.myVelocity(8,:),'r','LineWidth',1.5);
plot(engine.myTimes(1:end),engine.myVelocity(9,:),'k','LineWidth',1.5); 
plot(engine.myTimes(1:end),engine.myVelocity(10,:),'b','LineWidth',1.5);
legend('Velx-2','Vely-2','Velz-2');
xlabel('time (seconds)');
ylabel('Velocity (m/s)');
title('Velocity of Body2 vs Time') ;
subplot(3,1,3);hold on;
plot(engine.myTimes(1:end),engine.myAcceleration(8,:),'r','LineWidth',1.5); 
plot(engine.myTimes(1:end),engine.myAcceleration(9,:),'k','LineWidth',1.5); 
plot(engine.myTimes(1:end),engine.myAcceleration(10,:),'b','LineWidth',1.5);
legend('Accelx-2','Accely-2','Accelz-2');
xlabel('time (seconds)');
ylabel('Acceleration (m/s^2)');
title('Acceleration of Body2 vs Time') ;
end