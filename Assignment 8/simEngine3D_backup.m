classdef simEngine3D
    %%
    properties(GetAccess = 'public',SetAccess = 'public')
        %Logistical
        myProblemID % stores ADM info
        myProblemInfo % stores MDL info
        myProblemType %Kinematics or Dynamics?
        myStepSize % Step size in sim
        
        myFinalTime % Final time of sim
        myNbodies %# of bodies in simulation
        myNKC %Number of equations of constraint needed 7*nbodies
        myNACE %Number of algebraic constraint equations
        %Number of constraint equations (7*nbodies), 6DOF + 1 Euler Parameter normalization
        %General Simulation information
        myTimes %
        myTimeStep %
        myCurrentTime
        myPosition % Position-values of every DOF for all timesteps
        myqi %Position-values at current timestep
        myqdi
        myVelocity % Velocity-Value of every DOF at each timestep
        myAcceleration % Acceleration-Value of every DOF at each timestep
        myKinCon %
        myBodies %
        myFunTimes %Functions from every provided constraint in
        %'file.mdl' evaluated at all the times
        myJac % [p'*p-1;
        %                 obj.myJac(6*obj.myNbodies+i,7*i-3:7*i) = 2*p';
        myPhi % [Phi_K (nc); Phi_P (nb)]
        myNu %
        myGamma % [gamma cons (nc); gamma_P (nb)(-2pd'*pd];
        myFlags % 4x1
        
        %Dynamics / Inverse Dynamics
        myM
        myJ
        myLambda %Lagrange mult [lamb_p ;lamb] i.e.[p*lamb_P term in EOM; forces/torques]
        myFreact %7*nb X timeSteps
        
        
    end
    %%
    methods
        
        function obj = simEngine3D(inputACF,inputMDL)
            %Constructor for Kinematic Analysis.
            %Input: 'input.ACF' contains info on how to run the simulation
            % 'input.MDL' contains information on the bodies and constraints used in
            % simulation
            
            obj.myProblemID = loadjson(inputMDL);
            obj.myProblemInfo = loadjson(inputACF);
            obj.myProblemType = obj.myProblemInfo.simulation;
            obj.myTimeStep = 1;
            obj.myStepSize = obj.myProblemInfo.stepSize;
            obj.myFinalTime = obj.myProblemInfo.tend;
            obj.myTimes = 0:obj.myStepSize:...
                obj.myProblemInfo.outputSteps*obj.myStepSize;
            if strcmp(obj.myProblemType, 'Kinematics')
                obj.myCurrentTime = 0+ obj.myStepSize;
            else
                obj.myCurrentTime = 0;
                %                 obj.myForces = obj.myProblemID.forces;
            end
            obj.myNbodies = numel(obj.myProblemID.bodies);
            nb = obj.myNbodies;
            obj.myNKC = numel(obj.myProblemID.constraints);
            nc = obj.myNKC;
            obj.myNACE = obj.myNKC+obj.myNbodies;
            %Initialize position velocity and accel of bodies to 0.
            obj.myPosition = zeros(7*obj.myNbodies,...
                obj.myProblemInfo.outputSteps);
            ts = obj.myProblemInfo.outputSteps;
            obj.myVelocity = obj.myPosition;
            obj.myAcceleration = obj.myPosition;
            obj.myKinCon = obj.myProblemID.constraints;
            obj.myBodies = obj.myProblemID.bodies;
            obj.myPhi = zeros(obj.myNACE,1);
            obj.myJac = zeros(obj.myNACE,obj.myNACE);
            zer = zeros(obj.myNKC,numel(obj.myTimes));
            obj.myFunTimes = struct('funt',zer,'funDt',zer,'funDDt',zer);
            obj.myFlags = [1 1 0 0];
            obj = obj.getFunctions;
            obj.myM = zeros(obj.myNbodies,1);
            obj.myJ = zeros(3*obj.myNbodies,1);
            obj.myLambda = zeros(nc+nb,ts);
            obj.myFreact = zeros(7*nb,ts);
            
            for i = 1:obj.myNbodies
                %Creates structs for each body that contains a 3x3 Mass and
                %3x3 Moment of inertia matrix
                %                 Mi = obj.myBodies(i).mass;
                %                 J = obj.myBodies(i).jbar;
                %                 M = [Mi 0 0; 0 Mi 0; 0 0 Mi];
                %                 Jbar = [J(1) 0 0; 0 J(2) 0; 0 0 J(3)];
                %                 obj.myM{i} = struct('M', M);
                %                 obj.myJ{i} = struct('J', Jbar);
                obj.myM(i) = obj.myBodies(i).mass;
                obj.myJ(3*i-2:3*i) = obj.myBodies(i).jbar;
            end
            
        end
        %%
        function obj = runDynamics(obj)
            %6 Step plan: http://sbel.wisc.edu/Courses/ME751/2016/Documents/lecture1019.pdf
            %Slides 15 - 17
            maxit = 16;
            tol = 10^-3;
            h = obj.myStepSize;
            
            %Step -1: Solve for initial accelerations and lagrange multipliers
            
            nb = obj.myNbodies;
            for i = 1:nb
                obj.myPosition(7*i-6:7*i,1) = obj.myBodies(i).q0;
                obj.myVelocity(7*i-6:7*i,1) = obj.myBodies(i).qd0;
            end
            nc = numel(obj.myKinCon);
            bigM = diag(reshape(repmat(obj.myM,3,1),3*obj.myNbodies,1));
            bigP = zeros(nb,4*nb);
            bigJ = zeros(4*nb,4*nb);
            gamP = zeros(nb,1);
            gam = zeros(nc*nb,1);
            timestep = 1;
            obj = obj.getGammaJac(timestep);
            jac = obj.myJac;
            phiR = zeros(nc,3*nb);
            phiP = zeros(nc,4*nb);
            %This assembles all the submatrices used in building the linear
            %system.
            for i = 1:nb
                
                p = obj.myPosition(7*i-3:7*i,1);
                pd = obj.myVelocity(7*i-3:7*i,1);
                e = p(2:end);

                eTIL = [ 0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];

                G = [-e, -eTIL+p(1)*eye(3)]; 
                Jbar = diag(obj.myJ(3*i-2:3*i));
                bigJ(4*i-3:4*i,4*i-3:4*i) = 4*G'*Jbar*G;
                bigP(i,4*i-3:4*i) = p';
                
                %[Generalized body forces; gen body torques];
                %Only body forces added right now
                F = zeros(7*nb,1);
                grav = obj.myProblemID.gravity;
                F(7*i-6:7*i-4,1) = obj.myBodies(i).mass*grav;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Still need to implement Gen forces and
                %torques: SEE SLIDE 26 
                %http://sbel.wisc.edu/Courses/ME751/2016/Documents/lecture1005.pdf
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %F(4-7) = 8p'GdotTJbarGdotp
                F(7*i-3:7*i) = 8*pd'*G'*Jbar*G*pd;

                phiR(1:end,3*i-2:3*i) = jac(1:end-nb,7*i-6:7*i-4);
                phiP(1:end,4*i-3:4*i) = jac(1:end-nb,7*i-3:7*i);
            end
            gamP = obj.myGamma(nc+1:end);%REMOVED A -1*
            gam = obj.myGamma(1:nc);
            jacr = zeros(nc,3*nb);
            jacp = zeros(nc,4*nb);
            for k = 1:nb
               jacr(1:end,3*k-2:3*k) = jac(1:end-nb,7*k-6:7*k-4);
               jacp(1:end,4*k-3:4*k) = jac(1:end-nb,7*k-3:7*k);
            end
            bigLS = [ bigM, zeros(3*nb,5*nb),jacr';...
                zeros(4*nb,3*nb), bigJ, bigP',jacp';...
                zeros(nb,3*nb), bigP, zeros(nb,nb), zeros(nb,nc);...
                jacr,jacp, zeros(nc,nb+nc) ];
            b = [F; gamP'; gam'];
            x = bigLS\b;
            obj.myLambda(:,timestep) = x(end-nb-nc+1:end);
            lamb = obj.myLambda(nb+1:end,timestep);
            obj.myFreact(1:3*nb,timestep) = -phiR'*lamb;
            obj.myFreact(3*nb+1:end,timestep) = -phiP'*lamb;
            for k=1:nb
            obj.myAcceleration(7*k-6:7*k,timestep) = x(7*k-6:7*k);
            end
            
            obj.myTimeStep = 2;
            %Step 0:Prime it=0 step with previous values of accel and lambda
            for t = 2:numel(obj.myTimes)
                for it = 1:maxit
                    if it == maxit
                        disp('Maxout in Quasi-Newton iteration')
                        
                    end
                    %Step 1: Compute position and velocity using BDF and
                    %accelerations
                    if t==2 %BDF1 primer for BDF2
                        betah = h^2;
                        
                        if it == 1
                            qdd = obj.myAcceleration(:,t-1);
                            lamb = obj.myLambda(nb+1:end,t-1);
                            lambP = obj.myLambda(1:nb,t-1);
                        else
                            qdd = obj.myAcceleration(:,t);
                            lamb = obj.myLambda(nb+1:end,t);
                            lambP = obj.myLambda(1:nb,t);
                        end
                        %%%%%%%%%%
                        rdd = qdd(7*(1:nb)-6:7*(1:nb)-4);
                        pdd = qdd(7*(1:nb)-3:7*(1:nb));
                        r = obj.myPosition(7*(1:nb)-6:7*(1:nb)-4,1)+h^2*rdd;
                        p = obj.myPosition(7*(1:nb)-3:7*(1:nb),1)+h^2*pdd;
                        rd = obj.myVelocity(7*(1:nb)-6:7*(1:nb)-4,1)+h*rdd;
                        pd = obj.myVelocity(7*(1:nb)-3:7*(1:nb),1)+h*pdd;
                        %%%%%%%%%%%%
                    else %BDF2
                        %Calculate Position and Velocity
                        betah = ((2/3)*h)^2;
                        
                        if it == 1
                            qdd = obj.myAcceleration(:,t-1);
                            lamb = obj.myLambda(nb+1:end,t-1);
                            lambP = obj.myLambda(1:nb,t-1);
                        else
                            qdd = obj.myAcceleration(:,t);
                            lamb = obj.myLambda(nb+1:end,t);
                            lambP = obj.myLambda(1:nb,t);
                        end
                        %%%%%%%%%%%%%%%
                        rdd = qdd(7*(1:nb)-6:7*(1:nb)-4);
                        pdd = qdd(7*(1:nb)-3:7*(1:nb));
                        
                        Cr = (4/3)*obj.myPosition(7*(1:nb)-6:...
                            7*(1:nb)-4,t-1)-(1/3)*...
                            obj.myPosition(7*(1:nb)-6:7*(1:nb)-4,t-2);
                        Cp = (4/3)*obj.myPosition(7*(1:nb)-3:...
                            7*(1:nb),t-1)-(1/3)*...
                            obj.myPosition(7*(1:nb)-3:7*(1:nb),t-2);
                        Cdr =(4/3)*obj.myVelocity(7*(1:nb)-6:...
                            7*(1:nb)-4,t-1)-(1/3)*...
                            obj.myVelocity(7*(1:nb)-6:7*(1:nb)-4,t-2);
                        Cdp =(4/3)*obj.myVelocity(7*(1:nb)-3:...
                            7*(1:nb),t-1)-(1/3)*...
                            obj.myVelocity(7*(1:nb)-3:7*(1:nb),t-2);
                        r = Cr+((2/3).^2*h.^2)*rdd;
                        p = Cp+((2/3).^2*h.^2)*pdd;
                        rd = Cdr+h*rdd;
                        pd = Cdp+h*pdd;
                        %%%%%%%%%%%%%%%%%%
                    end
                    
                    %Save the calculated positions and accelerations.
                    obj.myPosition(7*(1:nb)-6:7*(1:nb)-4,t) = r;
                    obj.myPosition(7*(1:nb)-3:7*(1:nb),t) = p;
                    obj.myVelocity(7*(1:nb)-6:7*(1:nb)-4,t) = rd;
                    obj.myVelocity(7*(1:nb)-3:7*(1:nb),t) = pd;
                    obj = obj.getGammaJac(t);
                    jac = obj.myJac;
                    %This assembles all the submatrices used in building the linear
                    %system.
                    for i = 1:nb
                        
                        p = obj.myPosition(7*i-3:7*i,t);
                        e = p(2:end);
                        eTIL = [ 0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];
                        G = [-e, -eTIL+p(1)*eye(3)]; %%%%%%%Added - to G
                        Jbar = diag(obj.myJ(3*i-2:3*i));
                        bigJ(4*i-3:4*i,4*i-3:4*i) = 4*G'*Jbar*G;
                        bigP(i,4*i-3:4*i) = p';
                        
                        %[Generalized body forces; gen body torques];
                        %Only body forces added right now
                        F = zeros(7*nb,1);
                        grav = obj.myProblemID.gravity;
                        F(7*i-6:7*i-4,1) = obj.myBodies(i).mass*grav;
                        F(7*i-3:7*i) = 8*pd'*G'*Jbar*G*pd;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %Still need to implement Gen forces and
                        %torques
                        %INDEXES HERE ARE WRONG FOR MULTIBODY I THINK
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        gamP = obj.myGamma(nc+1:end);%REMOVED A -1*
                        gam = obj.myGamma(1:nc);
%                         gamP(i) = -1*obj.myGamma(7*i);
%                         gam(6*i-5:6*i) = obj.myGamma(7*i-6:7*i-1);
                        phiR(1:end,3*i-2:3*i) = jac(1:end-nb,7*i-6:7*i-4);
                        phiP(1:end,4*i-3:4*i) = jac(1:end-nb,7*i-3:7*i);
                    end
                    %Step 2: Compute residual of system g(z)
                    PSY = [ bigM, zeros(3*nb,5*nb),jac(1:end-nb,1:3)';...
                        zeros(4*nb,3*nb), bigJ, bigP',jac(1:end-nb,4:7)';...
                        zeros(nb,3*nb), bigP, zeros(nb,nb), zeros(nb,nc);...
                        jac(1:end-nb,:), zeros(nc,nb+nc) ];
                    z = [rdd; pdd; lambP; lamb];
                    gtop = [ bigM, zeros(3*nb,5*nb),jac(1:end-nb,1:3)';...
                        zeros(4*nb,3*nb), bigJ, bigP',jac(1:end-nb,4:7)']*z;
                    obj = obj.getPhiJac(t);
                    %myPhi is [Phi;Phi_pnorm], change order here to be in
                    %same form as procedure;
                    PHI = (1/betah)*[obj.myPhi(obj.myNKC+1:end);...
                        obj.myPhi(1:obj.myNKC)];
                    res = [gtop-F; PHI];
                    %Step 3: Solve BigPHI*delz=-g(z) to get delz.
                    %Step 4: z(v+1) = z(v) + delz(v)
                    %Step 5: v = v+1
                    %Step 6: If converged, recompute position and velocity and
                    %store. t = t+1, restart step 0.
                    delz = PSY\-res;
                    z = z+delz;
                    
                    obj.myAcceleration(7*(1:3)-6:7*(1:3)-4,t) = z(1:3*nb);
                    obj.myAcceleration(7*(1:3)-3:7*(1:3),t) = z(3*nb+1:7*nb);
                    obj.myLambda(:,t) = z(7*nb+1:end);
                    obj.myFreact(1:3*nb,t) = -phiR'*lamb;
                    obj.myFreact(3*nb+1:end,t) = -phiP'*lamb;
%                     cond(PSY)
%                     norm(delz)
                    
                    
                    
                    %Check if converged
                    if norm(delz)<=tol
                        if t == 2
                            r = obj.myPosition(7*(1:nb)-6:7*(1:nb)-4,1)+h^2*rdd;
                            p = obj.myPosition(7*(1:nb)-3:7*(1:nb),1)+h^2*pdd;
                            rd = obj.myVelocity(7*(1:nb)-6:7*(1:nb)-4,1)+h*rdd;
                            pd = obj.myVelocity(7*(1:nb)-3:7*(1:nb),1)+h*pdd;
                        else
                            r = Cr+((2/3).^2*h.^2)*z(1:3*nb);
                            p = Cp+((2/3).^2*h.^2)*z(3*nb+1:7*nb);
                            rd = Cdr+h*z(1:3*nb);
                            pd = Cdp+h*z(3*nb+1:7*nb);
                        end
                        %Save all variables
                        obj.myPosition(7*(1:3)-6:7*(1:3)-4,t) = r;
                        obj.myPosition(7*(1:3)-3:7*(1:3),t) = p;
                        
                        obj.myVelocity(7*(1:3)-6:7*(1:3)-4,t) = rd;
                        obj.myVelocity(7*(1:3)-3:7*(1:3),t) = pd;
                        break
                    end
                    %                         %Save all variables
                    %                         obj.myPosition(7*(1:3)-6:7*(1:3)-4,t) = r;
                    %                         obj.myPosition(7*(1:3)-3:7*(1:3),t) = p;
                    %
                    %                         obj.myVelocity(7*(1:3)-6:7*(1:3)-4,t) = rd;
                    %                         obj.myVelocity(7*(1:3)-3:7*(1:3),t) = pd;
                    %
                    %                         obj.myAcceleration(7*(1:3)-6:7*(1:3)-4,t) = z(1:3*nb);
                    %                         obj.myAcceleration(7*(1:3)-3:7*(1:3),t) = z(3*nb+1:7*nb);
                    %
                    %                         obj.myLambda = z(7*nb+1:end);
                end
                obj.myTimeStep = obj.myTimeStep+1;
            end
        end
        %%
        function obj = getFunctions(obj)
            %Really bad implementation to numerically estimate the time
            %derivatives of the Kin Driving constraints of the system.
            for i = 1:obj.myNKC
                h = 10^-4;
                delh = obj.myStepSize*h;
                t = 0:delh:(obj.myFinalTime+obj.myStepSize);
                funt = str2func(obj.myKinCon{i}.fun);
                funt = funt(t);
                funDt = diff(funt)/delh;
                funDDt = diff(funDt)/delh;
                obj.myFunTimes.funt(i,:) = funt(1:h^(-1):numel(funt(:))-h^(-1));
                obj.myFunTimes.funDt(i,:) = funDt(1:h^(-1):numel(funt(:))-h^(-1));
                obj.myFunTimes.funDDt(i,:) = funDDt(1:h^(-1):numel(funt(:))-h^(-1));
            end
        end
        %%
        function obj = PositionAnalysis(obj)
            % PositionAnalysis - This will populate the myPosition array
            % which includes the position information for all timesteps.
            % This is done through Newton Raphson iterative method for each
            % timestep.
            obj.myFlags = [1 1 0 0];
            tol = 10^-6;
            maxit = 10;
            
            %initial guess at time = 0;
            for j = 1:obj.myNbodies
                obj.myPosition(7*j-6:7*j,1) = obj.myBodies(j).q0;
            end
            q = obj.myPosition(:,1);
            for t = 1:obj.myFinalTime/obj.myStepSize
                for it = 1:maxit;
                    if it == maxit
                        disp('maxout')
                    end
                    obj = getPhiJac(obj, it);
                    corr = obj.myJac\obj.myPhi;
                    q = q - corr;
                    obj.myPosition(:,t) = q;
                    obj.myqi = q;
                    if norm(corr) <= tol
                        break;
                    end
                end
                obj.myTimeStep = obj.myTimeStep +1;
                obj.myCurrentTime = obj.myCurrentTime + obj.myStepSize;
            end
        end
        %%
        function obj = getPhiJac(obj,it)
            %Iterate through all given constraint equations (including
            %driving constraints given in the form of Kinematic
            %Constraints) and put them in the appropriate global matrices
            %myPhi and myJac
            
            obj.myFlags = [1,1,0,0];
            
            for i = 1:obj.myNKC
                b1 = obj.myKinCon{i}.body1;
                b2 = obj.myKinCon{i}.body2;
                fh = str2func(obj.myKinCon{i}.type);
                
                if it==1 && obj.myTimeStep~=1
                    timeStep = obj.myTimeStep-1;
                else
                    timeStep = obj.myTimeStep;
                end
                
                %if body 2 is ground, ignore the outputs from the
                %constraint functions
                if strcmp(b2,'ground')
                    obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        0;0;0;1;0;0;0];
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        0;0;0;0;0;0;0];
                    [PHI,JAC,~,~] = fh(obj.myKinCon{i},obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myPhi(i) = PHI;
                    obj.myJac(i,7*b1-6:7*b1) = JAC(1:7);
                elseif strcmp(b1,'ground')
                    obj.myqi = [0;0;0;1;0;0;0;...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [zeros(7,1);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [PHI,JAC,~,~] = fh(obj.myKinCon{i},obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myPhi(i) = PHI;
                    obj.myJac(i,7*b2-6:7*b2) = JAC(8:14);
                else
                    obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [PHI,JAC,~,~] = fh(obj.myKinCon{i},...
                        obj.myTimeStep,obj.myFunTimes,...
                        obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myPhi(i) = PHI;
                    obj.myJac(i,7*b1-6:7*b1) = JAC(1:7);
                    obj.myJac(i,7*b2-6:7*b2) = JAC(8:14);
                end
            end
            %Append all Euler parameter normalization constrants to the
            %end of obj.myPhi and obj.myJac
            
            for i = 1:obj.myNbodies
                p = obj.myPosition(7*i-3:7*i,timeStep);
                obj.myPhi(obj.myNKC+i) = (0.5)*(p'*p-1);
                obj.myJac(obj.myNKC+i,7*i-3:7*i) = 2*p';
            end
        end
        %%
        function obj = VelocityAnalysis(obj)
            obj.myNu = zeros(obj.myNACE,1);
            obj.myTimeStep = 1;
            obj.myCurrentTime = 0+ obj.myStepSize;
            obj.myFlags = [0,1,1,0];
            for j = 1:obj.myNbodies
                obj.myVelocity(7*j-6:7*j,1) = obj.myBodies(j).qd0;
            end
            
            for timestep = 1:obj.myFinalTime/obj.myStepSize
                obj = getNuJac(obj);
                obj.myVelocity(:,timestep) = obj.myJac\obj.myNu;
                obj.myTimeStep = obj.myTimeStep +1;
                obj.myCurrentTime = obj.myCurrentTime + obj.myStepSize;
            end
        end
        %%
        function obj = getNuJac(obj)
            obj.myFlags = [0,1,1,0];
            
            for i = 1:obj.myNKC
                b1 = obj.myKinCon{i}.body1;
                b2 = obj.myKinCon{i}.body2;
                fh = str2func(obj.myKinCon{i}.type);
                timeStep = obj.myTimeStep;
                %if body 2 is ground, ignore the outputs from the
                %constraint functions
                if strcmp(b2,'ground')
                    obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        0;0;0;1;0;0;0];
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        0;0;0;0;0;0;0];
                    [~,JAC,NU,~] = fh(obj.myKinCon{i},obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myNu(i) = NU;
                    obj.myJac(i,7*b1-6:7*b1) = JAC(1:7);
                elseif strcmp(b1,'ground')
                    obj.myqi = [0;0;0;1;0;0;0;...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [zeros(7,1);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [~,JAC,NU,~] = fh(obj.myKinCon{i},obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myNu(i) = NU;
                    obj.myJac(i,7*b2-6:7*b2) = JAC(8:14);
                    
                else
                    obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [~,JAC,NU,~] = fh(obj.myKinCon{i},obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myNu(i) = NU;
                    obj.myJac(i,7*b1-6:7*b1) = JAC(1:7);
                    obj.myJac(i,7*b2-6:7*b2) = JAC(8:14);
                end
            end
            %Append all Euler parameter normalization constrants to the
            %end of obj.myNu and obj.myJac
            for i = 1:obj.myNbodies
                p = obj.myPosition(7*i-3:7*i,timeStep);
                obj.myNu(6*obj.myNbodies+i) = 0;
                obj.myJac(6*obj.myNbodies+i,7*i-3:7*i) = 2*p';
            end
        end
        %%
        function obj = AccelerationAnalysis(obj)
            obj.myGamma = zeros(obj.myNKC+obj.myNbodies,1);
            obj.myTimeStep = 1;
            obj.myCurrentTime = 0 + obj.myStepSize;
            obj.myFlags = [0,1,0,1];
            
            for timestep = 1:obj.myFinalTime/obj.myStepSize
                obj = obj.getGammaJac();
                obj.myAcceleration(:,timestep) = obj.myJac\obj.myGamma;
                obj.myTimeStep = obj.myTimeStep +1;
                obj.myCurrentTime = obj.myCurrentTime + obj.myStepSize;
            end
        end
        %%
        function obj = getGammaJac(obj,it)
            obj.myFlags = [0,1,0,1];
            
            for i = 1:obj.myNKC
                b1 = obj.myKinCon{i}.body1;
                b2 = obj.myKinCon{i}.body2;
                fh = str2func(obj.myKinCon{i}.type);
                if (nargin==2)
                    timeStep = it;
                else
                    timeStep = obj.myTimeStep;
                end
                
                %if body 2 is ground, ignore the outputs from the
                %constraint functions
                if strcmp(b2,'ground')
                    obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        0;0;0;1;0;0;0];
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        0;0;0;0;0;0;0];
                    [~,JAC,~,GAMMA] = fh(obj.myKinCon{i},obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myGamma(i) = GAMMA;
                    obj.myJac(i,7*b1-6:7*b1) = JAC(1:7);
                elseif strcmp(b1,'ground')
                    obj.myqi = [0;0;0;1;0;0;0;...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [zeros(7,1);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [~,JAC,~,GAMMA] = fh(obj.myKinCon{i},obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myGamma(i) = GAMMA;
                    obj.myJac(i,7*b2-6:7*b2) = JAC(8:14);

                else
                    obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [~,JAC,~,GAMMA] = fh(obj.myKinCon{i},obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myGamma(i) = GAMMA;
                    obj.myJac(i,7*b1-6:7*b1) = JAC(1:7);
                    obj.myJac(i,7*b2-6:7*b2) = JAC(8:14);
                end
            end
            %Append all Euler parameter normalization constrants to the
            %end of obj.myGamma and obj.myJac
            for i = 1:obj.myNbodies
                p = obj.myPosition(7*i-3:7*i,timeStep);
                pdi = obj.myVelocity(7*i-3:7*i,timeStep);
                obj.myGamma(obj.myNKC+i) = -2*(pdi'*pdi);
                obj.myJac(obj.myNKC+i,7*i-3:7*i) = 2*p';
            end
        end
        %%
        function obj = InverseDynamics(obj)
            obj.myTimeStep = 1;
            nb = obj.myNbodies;
            nc = obj.myNbodies*6;
            bigM = diag(reshape(repmat(obj.myM,3,1),3*obj.myNbodies,1));
            bigP = zeros(nb,4*nb);
            bigJ = zeros(4*nb,4*nb);
            gamP = zeros(nb,1);
            gam = zeros(nc*nb,1);
            %'the big linear system', Slide 4 of:
            %http://sbel.wisc.edu/Courses/ME751/2016/Documents/lecture1007.pdf
            
            for timestep = 1:obj.myFinalTime/obj.myStepSize
                obj = obj.getGammaJac(timestep);
                jac = obj.myJac;
                phiR = zeros(nc,3*nb);
                phiP = zeros(nc,4*nb);
              
                for i = 1:nb
                    p = obj.myPosition(7*i-3:7*i);
                    e = p(2:end);
                    eTIL = [ 0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];
                    G = [-e', -eTIL+p(1)*eye(3)];
                    Jbar = diag(obj.myJ(3*i-2:3*i));
                    bigJ(4*i-3:4*i,4*i-3:4*i) = 4*G'*Jbar*G;
                    bigP(i,4*i-3:4*i) = p';
                    %[Generalized body forces; gen body torques];
                    %Still need to be implemented
                    F = zeros(7*nb,1);
                    % INDEX IS WRONG HERE FOR MULTIBODY
                    gamP(i) = -1*obj.myGamma(obj.myNKC+1:end);
                    gam(6*i-5:6*i) = obj.myGamma(7*i-6:7*i-1);
                    phiR(1:end,3*i-2:3*i) = jac(1:end-1,7*i-6:7*i-4);
                    phiP(1:end,4*i-3:4*i) = jac(1:end-1,7*i-3:7*i);
 
                end
                bigLS = [ bigM, zeros(3*nb,5*nb),jac(1:end-1,1:3)';...
                    zeros(4*nb,3*nb), bigJ, bigP',jac(1:end-1,4:7)';...
                    zeros(nb,3*nb), bigP, zeros(nb,nb), zeros(nb,nc);...
                    jac(1:end-1,:), zeros(nc,nb+nc) ];
                b = [F; gamP; gam];
                cond(bigLS)
                x = bigLS\b;
                obj.myLambda(:,timestep) = x(end-nb-nc+1:end);
                lamb = obj.myLambda(nb+1:end,timestep);
                obj.myFreact(1:3*nb,timestep) = -phiR'*lamb;
                obj.myFreact(3*nb+1:end,timestep) = -phiP'*lamb;
            end
        end
        %%
        
        
    end
    
end