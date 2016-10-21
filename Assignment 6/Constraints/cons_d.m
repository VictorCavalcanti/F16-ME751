function [ Phi, Phi_q, Nu, Gamma] = cons_d(constraint, time, qi, qdi, flags)
%% -------------------------------------------------------------------------
% INPUTS:
% -constraint- [struct] Built from "bodies_info.adm".
% {fields} -  {"id": 4,"type": "D","body1": 1,
%               "sPiBAR": [ai_x, ai_y, ai_z],"body2": "2",
%               "sQjBAR":[aj_x, aj_y, aj_z],"fun": "f(t)^2"}
% -time- [scalar] current time
% -qi- [14x1 vector] containing: 
% [r1i, r2i, r3i, p0i, p0i, p1i, p2i, p3i, r1j, r2j, r3j, p0j, p1j, p2j,
%  p3j]
% -qdi- [14x1 vector] containing the time derivative of the previous
% quantities.
% -flags- [4x1 vector] controls what the function returns.
% 1-return; 0-no return; [Phi, Phi_q, Nu, Gamma]
% OUTPUTS:
% Phi - [Scalar]
% Phi_q - [1x14 vector] Jacobian of Phi w.r.t qi then qj
% Nu - -f_dot(t)
% Gamma - [scalar] Gamma = -(Phi_q*qdi)_q*qdi-2*Phi_q*qdi-Phi_tt
%% -------------------------------------------------------------------------
ri = qi(1:3); rid = qdi(1:3);
pi = qi(4:7); pid = qdi(4:7);
ei = pi(2:4); 
% eid = qdi(2:4);
eiTIL = [ 0 -ei(3) ei(2); ei(3) 0 -ei(1); -ei(2) ei(1) 0];
Ai = (pi(1).^2-(ei')*ei)*eye(3)+2*ei*(ei')+2*pi(1)*eiTIL;
sPiBAR = constraint.sPiBAR';
rj = qi(7+(1:3)); rjd = qdi(8:10);
pj = qi(7+(4:7)); pjd = qdi(11:14);
ej = pj(2:4); 
% ejd = qdj(2:4);
ejTIL = [ 0 -ej(3) ej(2); ej(3) 0 -ej(1); -ej(2) ej(1) 0];
Aj = (pj(1).^2-(ej')*ei)*eye(3)+2*ej*(ei')+2*pj(1)*ejTIL;
sQjBAR = constraint.sQjBAR';
dij = rj+Aj*sQjBAR-ri-Ai*sPiBAR;

%% Initialize outputs to empty, in case they will not be calculated.
Phi = [];
Phi_q = [];
Nu = [];
Gamma = [];
%% Calculate time derivatives used in this constraint:
t = sym('t');
fun_sym = sym(constraint.fun);
funt = matlabFunction(fun_sym, 'vars', t);
funDt = matlabFunction(diff(fun_sym), 'vars', t);
funDDt = matlabFunction(diff(diff(fun_sym)),'vars',t);


%% Calculate Phi - [scalar]
if flags(1)
    Phi = dij'*dij - funt(time);
end
%% Calculate Jacobian of Phi-[1x14 vector]
if flags(2)
    B_pi_spi = getB(pi,sPiBAR);
    B_pj_sqj = getB(pj,sQjBAR);
    
    %ADD A CHECK FOR GROUND HERE? 
    Phi_q = [-2*dij', -2*dij'*B_pi_spi, 2*dij', 2*dij'*B_pj_sqj];
%     Phi_q = [-c', -c'*B_pi_sPi, c', c'*B_pj_sjb];
    
end
%% Calculate Nu - [Scalar] Found from: Phi_q*qddi = -Phi_t = Nu    
if flags(3)
    Nu = funDt(time);
end
%% Calculate Gamma - [scalar]
if flags(4)
    B_pid_spi = getB(pid,sQjBAR);
    B_pjd_sqj = getB(pjd,sQjBAR);
    dijd = rjd + B_pj_sqj*pjd - rid - B_pi_spi*pid;
    Gamma = -2*dij'*B_pjd_sqj*pjd + 2*dij'*B_pid_spi*pid - 2*(dijd'*dijd)...
            +funDDt(time);
end
%% Calculate B
    function [B] = getB(p,a)
        e = p(2:4); 
        e0 = p(1);
        aTIL = [ 0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
        eTIL = [ 0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];
        B = [2*((e0.^2)*eye(3)+eTIL)*a, ...
            2*(e*a'-(e0*eye(3)+eTIL)*aTIL)];
    end
end
