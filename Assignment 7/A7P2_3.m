% IVP (RHS + IC)
clear all; clc; close all;
f = @(lam,y) lam*y;
y0 = 1;
tend = .02;
lam = -1000;
% Analytical solution
figure, hold on, box on
h = [0.001 0.0015 0.002];
colors = [[0, 0.4, 0]; [1, 0.5, 0]; [0.6, 0, 0]];
for ih = 1:length(h)
    tspan = 0:h(ih):tend;
    y_an = exp(lam'*tspan);
    y = zeros(3,numel(tspan));
    %     err = zeros(size(tspan));
    y(:,1) = y0;
    %     err(1) = 0;
    for i = 2:length(tspan)
        y(ih,i) = y(ih,i-1) + h(ih) * f(lam, y(ih,i-1));
        % err(i) = y(i) - y_an(tspan(i));
        
    end
    plot(tspan, y(ih,:), 'color', colors(ih,:),'lineWidth', 1.5);
    
    
end
plot(tspan,y_an(:));
legend('h = 0.001', 'h = 0.0015', 'h = 0.002');
title('Lambda = -1000');