function [ Phi, Phi_q, Nu, Gamma] = cons_cd(constraint, time,funtimes, qi, qdi, flags)
%% -------------------------------------------------------------------------
% INPUTS:
% -constraint- [struct] Built from "bodies_info.adm".
% {fields} -  {"id": 1,"type": "CD","c": [1, 0, 0],"body1": 1,
%               "sPi": [ai_x, ai_y, ai_z],"body2": "2", "sQj":
%               [aj_x, aj_y, aj_z],"fun": "f(t)"}
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
c = constraint.c';
ri = qi(1:3);
pi = qi(4:7); pid = qdi(4:7);
ei = pi(2:4);
eiTIL = [ 0 -ei(3) ei(2); ei(3) 0 -ei(1); -ei(2) ei(1) 0];
Ai = (pi(1).^2-(ei')*ei)*eye(3)+2*ei*(ei')+2*pi(1)*eiTIL;
sPiBAR = constraint.sPiBAR';
rj = qi(7+(1:3));
pj = qi(7+(4:7)); pjd = qdi(11:14);
ej = pj(2:4);
ejTIL = [ 0 -ej(3) ej(2); ej(3) 0 -ej(1); -ej(2) ej(1) 0];
Aj = (pj(1).^2-(ej')*ej)*eye(3)+2*ej*(ej')+2*pj(1)*ejTIL;
sQjBAR = constraint.sQjBAR';
id = constraint.id;
%% Initialize outputs to empty, in case they will not be calculated.
Phi = [];
Phi_q = [];
Nu = [];
Gamma = [];

%% Calculate Phi - [scalar]
if flags(1)
    Phi = c'*(rj+Aj*sQjBAR-ri-Ai*sPiBAR) - funtimes.funt(id,time);
end
%% Calculate Jacobian of Phi-[1x14 vector]
if flags(2)
    B_pi_sPi = getB(pi,sPiBAR);
    B_pj_sQj = getB(pj,sQjBAR);
    Phi_q = [-c', -c'*B_pi_sPi, c', c'*B_pj_sQj];
    
end
%% Calculate Nu - [Scalar] Found from: Phi_q*qddi = -Phi_t = Nu
if flags(3)
Nu = +funtimes.funDt(id,time);
end
%% Calculate Gamma - [scalar]
if flags(4)
    B_pjd_sQj = getB(pjd,sQjBAR);
    B_pid_sPi = getB(pid,sPiBAR);
    Gamma = -c'*B_pjd_sQj*pjd + c'*B_pid_sPi*pid + ...
        funtimes.funDDt(id,time);
end
%% Calculate B
    function [B] = getB(p,a)
        e = p(2:4);
        e0 = p(1);
        aTIL = [ 0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
        eTIL = [ 0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];
        B = [2*(e0*eye(3)+eTIL)*a, ...
            2*(e*a'-(e0*eye(3)+eTIL)*aTIL)];
    end
end
