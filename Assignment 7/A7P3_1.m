close all; clear all; clc;
for i = 1:4
alph = [0 1 2 3];
beta = 1;
x0 = 0;
y0 = 2;
tend = 20;
h = 0.05;
tspan = 0:h:tend;
x = zeros(1,length(tspan));
y = zeros(1,length(tspan));
y(1) = 2;
y(2) = 2;
tol = 10^-6;
maxit = 10;
% colors = [[1, 0, 0]; [0, 0, 0];[ 0, 0, 1];  [0, 1, 0]];
colors = [ 'r--'; 'g--'; 'k--'; 'b--'];
colors2 = [ 'r-'; 'g-'; 'k-'; 'b-'];
for t = 1:length(tspan)-1
    q = [x(t); y(t)];
    for it = 1:maxit;
        if it == maxit
            disp('maxout')
        end
        g = [x(t+1)*(1+h)+((4*h*x(t+1)*y(t+1))/(1+x(t+1)^2))-x(t)-h*alph(i);...
            -h*beta*x(t+1)+y(t+1)+(h*beta*x(t+1)*y(t+1))/(1+x(t+1)^2)-y(t)];
        Jac = [ 1+h+4*h*y(t+1)*(1-x(t+1)^2)/(1+x(t+1)^2).^2, 4*h*x(t+1)/(1+x(t+1)^2);...
        -h*beta+h*beta*(1-x(t+1)^2)/(1+x(t+1)^2)^2 , 1+(beta*h*x(t+1))/(1+x(t+1)^2)];
        corr = Jac\(g);
        q = q - corr;
        x(t+1) = q(1);
        y(t+1) = q(2);
        if norm(corr) <= tol
            break;
        end
    end
end
plot(tspan,x,colors(i,:),'lineWidth',1.5); hold on;
plot(tspan,y,colors2(i,:),'lineWidth',1.5);
legend('alpha = 0, x','alpha = 0, y','alpha = 1, x','alpha = 1, y',...
    'alpha = 2, x','alpha = 2, y','alpha = 3, x','alpha = 3, y' );
title('Different values of alpha for beta = 1');

end