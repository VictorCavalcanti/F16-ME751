close all; clear all; clc;
tend = 10;
y1 = 1;
tol = 10^-6;
hLargest = .01;
maxit = 15;
nPoints = 8;
hSize = zeros(8,1);
hSize(1) = hLargest;

for i=1:nPoints-1
    hSize(i+1)=hSize(i)/2;
end
results = zeros(nPoints, 3);
for j = 1:nPoints
    tspan = 1:hSize(j):tend;
    y_an = (1./tspan)+(1./tspan.^2).*tan((1./tspan)+pi-1);
    y = zeros(size(tspan));
    y2 = zeros(size(tspan));
    y(1) = y1;
    y(2) = y1;
    y2(1:5) = y_an(1:5);
    for t = 1:length(tspan)-1
        q = y(t);
        q2 = y2(t);
        for it = 1:maxit;
            if it == maxit
                disp('maxout')
            end
            g = y(t+1)+hSize(j)*y(t+1)^2+(hSize(j)/tspan(t+1)^4)-y(t);
            %             Jac = 1+2*hSize(j)*y(t+1)-4*hSize(j)/tspan(t+1)^5;
            
            Jac = 1+2*hSize(j)*y(t+1);
            corr = Jac\(g);
            q = q - corr;
            y(t+1) = q;
            if norm(corr) <= tol
                break;
            end
        end
        for it = 1:maxit;
            
            if t >= 4
                if it==maxit
                    disp('maxoutBDF')
                end
                gBDF = y2(t+1)-(48/25)*y2(t)+(36/25)*y2(t-1)-(16/25)*y2(t-2)+...
                    (3/25)*y2(t-3)-(12/25)*hSize(j)*(-y2(t+1).^2-(1/tspan(t+1).^4));
                Jac2 = 1+(24/25)*hSize(j)*y2(t+1);
                corr2 = Jac2\gBDF;
                q2 = q2 - corr2;
                y2(t+1) = q2;
                if norm(corr2) <= tol
                    break;
                end
            end
        end
    end
    results(j,1) = hSize(j);
    results(j,2) = abs(y_an(end)-y(end));
    results(j,3) = abs(y_an(end)-y2(end));
end
plot(log2(results(:,1)),log2(results(:,2)),'lineWidth',1.5); hold on;
plot(log2(results(:,1)),log2(results(:,3)),'k-','lineWidth',1.5);
legend('y_an(end)-y_BE(end)', 'y_an(end)-y_BDF(end)');
xlabel('log2(hSize)');
ylabel('log2(abs(y_soln-y_n))');
slope1 = (results(5,2)-results(4,2))/(results(5,1)-results(4,1));
slope2 = (results(5,3)-results(4,3))/(results(5,1)-results(4,1));