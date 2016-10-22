    function [B] = getB(p,a)
        e = p(2:4);
        e0 = p(1);
        aTIL = [ 0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
        eTIL = [ 0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];
        B = [2*(e0*eye(3)+eTIL)*a, ...
            2*(e*a'-(e0*eye(3)+eTIL)*aTIL)];
    end