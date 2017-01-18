function [A] = getA(p)
%found in Key Kinematic Equations:
% http://sbel.wisc.edu/Courses/ME751/2010/Documents/kinematicsKeyFormulas.pdf
A = [ p(1).^2+p(2).^2-p(3).^2-p(4).^2 , 2*(p(2)*p(3)-p(1)*p(4)),...
    2*(p(2)*p(4)+p(1)*p(3)); 2*(p(2)*p(3)+p(1)*p(4)) ...
p(1).^2-p(2).^2+p(3).^2-p(4).^2 , 2*(p(3)*p(4)-p(1)*p(2));...
2*(p(2)*p(4)-p(1)*p(3)) , 2*(p(3)*p(4)+p(1)*p(2)), p(1).^2-p(2).^2-p(3).^2+p(4).^2];
end
