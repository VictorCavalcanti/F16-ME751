classdef simEngine3D
    %%
    properties(GetAccess = 'public',SetAccess = 'public')
        myTemp %TEMP
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
        myTimes % vector with the actual time at east timestep
        myTimeStep % current timestep
        myCurrentTime % current time
        myPosition % Position-values of every DOF for all timesteps
        myqi %Position-values at current timestep
        myqdi
        myVelocity % Velocity-Value of every DOF at each timestep
        myAcceleration % Acceleration-Value of every DOF at each timestep
        myKinCon % Struct with constraint information
        myBodies % Struct with information about the bodies
        myFunTimes %Functions from every provided constraint in
        %'file.adm' evaluated at all the times
        myJac % [p'*p-1;
        myNumJac %Jacobian computed via finite difference.
        myPhi % [Phi_K (nc); Phi_P (nb)]
        myNu %
        myGamma % [gamma cons (nc); gamma_P (nb)(-2pd'*pd];
        myFlags % 5x1 flags for calculating quantities in the constraints
        myITS % 1xtimeSteps # iterations for convergence
        %Dynamics / Inverse Dynamics
        myLambda %Lagrange mult [lamb_p ;lamb] i.e.[p*lamb_P term in EOM; forces/torques]
        myFreact %7*nb X timeSteps
        myM
        myJ
        %Matrices used in calculating the Newton-Raphson Jacobian
        myBigM %3nbx3nb mass matrix
        myBigP %3nbx4nb Euler param matrix
        myBigJ
        myNRIjac %Struct containig fields (body1, body2, jacobian);
        myPSY %Jacobian used in Newton Raphson iterations.
        myBetaH %Term that arises from using BDF for NR
        myTSDA %Translational Spring Damper Actuator
        myNTSDA %Number of TSDAs
        myBDF6ON %Flag to see if solver should use BDF order 6.
        
        
        
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
            for i = 1:numel(obj.myKinCon)
               obj.myKinCon{i}.id = i; 
            end
            obj.myBodies = obj.myProblemID.bodies;
            obj.myPhi = zeros(obj.myNACE,1);
            obj.myJac = zeros(obj.myNACE,obj.myNACE);
            zer = zeros(obj.myNKC,numel(obj.myTimes));
            obj.myFunTimes = struct('funt',zer,'funDt',zer,'funDDt',zer);
            obj.myFlags = [1 1 0 0 0];
            obj = obj.getFunctions;
            obj.myM = zeros(obj.myNbodies,1);
            obj.myJ = zeros(3*obj.myNbodies,1);
            obj.myLambda = zeros(nc+nb,ts);
            obj.myFreact = zeros(7*nb,ts);
            obj.myITS = zeros(1,ts);
            obj.myTemp = zeros(nb,ts); % TEMP FOR HW8
            obj.myTSDA = obj.myProblemID.actuators;
            obj.myNTSDA = numel(obj.myTSDA);
            for i = 1:obj.myNbodies
                %Creates structs for each body that contains a 3x3 Mass and
                %3x3 Moment of inertia matrix
                obj.myM(i) = obj.myBodies(i).mass;
                obj.myJ(3*i-2:3*i) = obj.myBodies(i).jbar;
            end
            
        end
        %%
        function obj = runDynamics(obj)
            %6 Step plan: http://sbel.wisc.edu/Courses/ME751/2016/Documents/lecture1019.pdf
            %Slides 15 - 17
            maxit = 50;
            tol = 10^-4;
            h = obj.myStepSize;
            
            %Step -1: Solve for initial accelerations and lagrange multipliers
            
            nb = obj.myNbodies;
            rindx = zeros(3*nb,1);
            pindx = zeros(4*nb,1);
            for i = 1:nb
                obj.myPosition(7*i-6:7*i,1) = obj.myBodies(i).q0;
                obj.myVelocity(7*i-6:7*i,1) = obj.myBodies(i).qd0;
                rindx(3*i-2:3*i) = (7*i-6:7*i-4);
                pindx(4*i-3:4*i) = (7*i-3:7*i);
            end
            nc = numel(obj.myKinCon);
            bigM = diag(reshape(repmat(obj.myM,1,3)',3*obj.myNbodies,1));
            bigP = zeros(nb,4*nb);
            bigJ = zeros(4*nb,4*nb);
            
            timestep = 1;
            obj = obj.getGammaJac(timestep);
            jac = obj.myJac;
            phiR = zeros(nc,3*nb);
            phiP = zeros(nc,4*nb);
            F = zeros(7*nb,1);
            
            %This assembles all the submatrices used in building the linear
            %system.
            for i = 1:nb
                p = obj.myPosition(7*i-3:7*i,1);
                pd = obj.myVelocity(7*i-3:7*i,1);
                G = obj.getG(p);
                Gd = obj.getG(pd);
                Jbar = diag(obj.myJ(3*i-2:3*i));
                bigJ(4*i-3:4*i,4*i-3:4*i) = 4*G'*Jbar*G;
                bigP(i,4*i-3:4*i) = p';
                %[Generalized body forces; gen body torques];
                %Only body forces added right now
                grav = obj.myProblemID.gravity;
                F(3*i-2:3*i) = obj.myBodies(i).mass*grav;
                F(3*nb+(4*i-3):(3*nb+(4*i))) = 8*Gd'*Jbar*Gd*p;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Still need to implement Gen forces and
                %torques: SEE SLIDE 26
                %http://sbel.wisc.edu/Courses/ME751/2016/Documents/lecture1005.pdf
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                phiR(1:end,3*i-2:3*i) = jac(1:end-nb,7*i-6:7*i-4);
                phiP(1:end,4*i-3:4*i) = jac(1:end-nb,7*i-3:7*i);
            end
            if obj.myNTSDA~=0
                ftsda = obj.getFTSDA();
                F = F+ftsda;
            end
            
            gamP = obj.myGamma(nc+1:end);
            gam = obj.myGamma(1:nc);
            jacr = jac(1:end-nb,rindx);
            jacp = jac(1:end-nb,pindx);
            bigLS = [ bigM, zeros(3*nb,5*nb),jacr';...
                zeros(4*nb,3*nb), bigJ, bigP',jacp';...
                zeros(nb,3*nb), bigP, zeros(nb,nb), zeros(nb,nc);...
                jacr, jacp, zeros(nc,nb+nc) ];
            b = [F; gamP'; gam'];
            x = bigLS\b;
            obj.myLambda(:,timestep) = x(end-nb-nc+1:end);
            lamb = obj.myLambda(nb+1:end,timestep);
            obj.myFreact(1:3*nb,timestep) = -phiR'*lamb;
            obj.myFreact(3*nb+1:end,timestep) = -phiP'*lamb;
            obj.myAcceleration([rindx; pindx],timestep) = x(1:7*nb);
%             pd = obj.myVelocity(pindx,1);
%             obj.myTemp(:,1) = bigP*pd; %P*pd=0nb lvl1 constraint
            %%
            obj.myTimeStep = 2;
            %Step 0:Prime it=0 step with previous values of accel and lambda
            for t = 2:numel(obj.myTimes)
                for it = 1:maxit
                    if it == maxit
                        disp('Maxout in Newton-Raphson iteration')
                    end
                    %Step 1: Compute position and velocity using BDF and
                    %accelerations
                    if it == 1
                        obj.myAcceleration(:,t) = ...
                            obj.myAcceleration(:,t-1);
                        obj.myLambda(:,t) = obj.myLambda(:,t-1);
                    end
                        qdd = obj.myAcceleration(:,t);
                        lamb = obj.myLambda(nb+1:end,t);
                        lambP = obj.myLambda(1:nb,t);
                        rdd = qdd(rindx);
                        pdd = qdd(pindx);
                    if t==2 %BDF1 primer for BDF2
                        obj.myBetaH = h^2;
                        if it == 1
                            obj.myAcceleration(:,t) = ...
                                obj.myAcceleration(:,t-1);
                            obj.myLambda(:,t) = obj.myLambda(:,t-1);
                        end

                        r = obj.myPosition(rindx)+h*obj.myVelocity(rindx)+h^2*rdd;
                        p = obj.myPosition(pindx)+h*obj.myVelocity(pindx)+h^2*pdd;
                        rd = obj.myVelocity(rindx)+h*rdd;
                        pd = obj.myVelocity(pindx)+h*pdd;
                        


                        
                    elseif t==3  %BDF2 primer for BDF3
                        %Calculate Position and Velocity
                        obj.myBetaH = ((2/3)*h)^2;

                        
                        Cr = (4/3)*obj.myPosition(rindx,t-1)-(1/3)*...
                            obj.myPosition(rindx,t-2)+...
                            (8/9)*h*obj.myVelocity(rindx,t-1)-...
                            (2/9)*h*obj.myVelocity(rindx,t-2);
                        Cp = (4/3)*obj.myPosition(pindx,t-1)-(1/3)*...
                            obj.myPosition(pindx,t-2)+...
                            (8/9)*h*obj.myVelocity(pindx,t-1)...
                            -(2/9)*h*obj.myVelocity(pindx,t-2);
                        Cdr =(4/3)*obj.myVelocity(rindx,t-1)-(1/3)*...
                            obj.myVelocity(rindx,t-2);
                        Cdp =(4/3)*obj.myVelocity(pindx,t-1)-(1/3)*...
                            obj.myVelocity(pindx,t-2);
                        r = Cr+obj.myBetaH*rdd;
                        p = Cp+obj.myBetaH*pdd;
                        rd = Cdr+sqrt(obj.myBetaH)*rdd;
                        pd = Cdp+sqrt(obj.myBetaH)*pdd;
                        

                      elseif strcmp(obj.myProblemInfo.BDF6ON,'yes') && t>7
                        obj.myBetaH = ((60/147)*h)^2;
                        rt = obj.myPosition(rindx,t-6:t-1);
                        pt = obj.myPosition(pindx,t-6:t-1);
                        vt = obj.myVelocity(rindx,t-6:t-1);
                        vpt = obj.myVelocity(pindx,t-6:t-1);
                        rt = fliplr(rt);
                        pt = fliplr(pt);
                        vt = fliplr(vt);
                        vpt = fliplr(vpt);

                        
                        Cr = (360/147)*rt(:,1)-(450/147)*rt(:,2)+(400/147)*rt(:,3)-(225/147)*rt(:,4)+(72/147)*rt(:,5)-(10/147)*rt(:,6)+h*(60/147)*((360/147)*vt(:,1)-(450/147)*vt(:,2)+(400/147)*vt(:,3)-(225/147)*vt(:,4)+(72/147)*vt(:,5)-(10/147)*vt(:,6));
                        Cp = (360/147)*pt(:,1)-(450/147)*pt(:,2)+(400/147)*pt(:,3)-(225/147)*pt(:,4)+(72/147)*pt(:,5)-(10/147)*pt(:,6)+h*(60/147)*((360/147)*vpt(:,1)-(450/147)*vpt(:,2)+(400/147)*vpt(:,3)-(225/147)*vpt(:,4)+(72/147)*vpt(:,5)-(10/147)*vpt(:,6));
                        Cdr = (360/147)*vt(:,1)-(450/147)*vt(:,2)+(400/147)*vt(:,3)-(225/147)*vt(:,4)+(72/147)*vt(:,5)-(10/147)*vt(:,6);
                        Cdp = (360/147)*vpt(:,1)-(450/147)*vpt(:,2)+(400/147)*vpt(:,3)-(225/147)*vpt(:,4)+(72/147)*vpt(:,5)-(10/147)*vpt(:,6);
                        r = Cr+obj.myBetaH*rdd;
                        p = Cp+obj.myBetaH*pdd;
                        rd = Cdr+sqrt(obj.myBetaH)*rdd;
                        pd = Cdp+sqrt(obj.myBetaH)*pdd;
                        
                      else %BDF3
                        %Calculate Position and Velocity
                        obj.myBetaH = ((6/11)*h)^2;

                        
                        Cr = (18/11)*obj.myPosition(rindx,t-1)-...
                            (9/11)*obj.myPosition(rindx,t-2)+...
                            (2/11)*obj.myPosition(rindx,t-3)+...
                            (6/11)*h*((18/11)*obj.myVelocity(rindx,t-1)-...
                            (9/11)*obj.myVelocity(rindx,t-2)+...
                            (2/11)*obj.myVelocity(rindx,t-3));
                        
                        Cp = (18/11)*obj.myPosition(pindx,t-1)-(9/11)*...
                            obj.myPosition(pindx,t-2)+...
                            (2/11)*obj.myPosition(pindx,t-3)+...
                            (6/11)*h*((18/11)*obj.myVelocity(pindx,t-1)-...
                            (9/11)*obj.myVelocity(pindx,t-2)+...
                            (2/11)*obj.myVelocity(pindx,t-3));

                        
                        Cdr = (18/11)*obj.myVelocity(rindx,t-1)-...
                            (9/11)*obj.myVelocity(rindx,t-2)+...
                            (2/11)*obj.myVelocity(rindx,t-3);
                        Cdp = (18/11)*obj.myVelocity(pindx,t-1)-...
                            (9/11)*obj.myVelocity(pindx,t-2)+...
                            (2/11)*obj.myVelocity(pindx,t-3);
                        r = Cr+obj.myBetaH*rdd;
                        p = Cp+obj.myBetaH*pdd;
                        rd = Cdr+sqrt(obj.myBetaH)*rdd;
                        pd = Cdp+sqrt(obj.myBetaH)*pdd;
                        
                    end

                    
                    
                    %Save the calculated positions and accelerations.
                    obj.myPosition([rindx; pindx],t) = [r ;p];
                    obj.myVelocity([rindx; pindx],t) = [rd; pd];
                    obj = obj.getGammaJac(t);
                    jac = obj.myJac;
                    jacr = jac(1:end-nb,rindx);
                    jacp = jac(1:end-nb,pindx);
                    %This assembles all the submatrices used in building the linear
                    %system.
                    F = zeros(7*nb,1);
                    %Calculate some forces and matrices used in linear
                    %system
                    for i = 1:nb
                        p = obj.myPosition(7*i-3:7*i,t);
                        pd = obj.myVelocity(7*i-3:7*i,t);
                        G = obj.getG(p);
                        Gd = obj.getG(pd);
                        Jbar = diag(obj.myJ(3*i-2:3*i));
                        bigJ(4*i-3:4*i,4*i-3:4*i) = 4*G'*Jbar*G;
                        bigP(i,4*i-3:4*i) = p';
                        %[Generalized body forces; gen body torques];
                        grav = obj.myProblemID.gravity;
                        F(3*i-2:3*i) = obj.myBodies(i).mass*grav;
                        F(3*nb+(4*i-3):(3*nb+(4*i))) = 8*Gd'*Jbar*Gd*p;
                        phiR(1:end,3*i-2:3*i) = jac(1:end-nb,7*i-6:7*i-4);
                        phiP(1:end,4*i-3:4*i) = jac(1:end-nb,7*i-3:7*i);
                    end
                    %ADD TSDA FORCES
                    if obj.myNTSDA~=0
                        ftsda = obj.getFTSDA();
                        F = F+ftsda;
                    end
                    obj.myBigM = bigM;
                    obj.myBigP = bigP;
                    obj.myBigJ = bigJ;
                    %%
                    %%%%Check to see which Newton Raphson Mothod to use:
                    %%%%Quasi (Simplified Jacobian being computed),
                    %%%%Modified (Full Jacobian, computed at it=1),
                    %%%%Full (Full Jacobian, computed at every it).
                    
                    switch obj.myProblemInfo.NRIT
                        case{'Quasi'}
                            obj.myPSY = [ bigM, zeros(3*nb,5*nb),jacr';...
                                zeros(4*nb,3*nb), bigJ, bigP',jacp';...
                                zeros(nb,3*nb), bigP, zeros(nb,nb), zeros(nb,nc);...
                                jacr, jacp, zeros(nc,nb+nc)];
                        case{'Modified'}
                            if it == 1
                                obj.myPSY = [ bigM, zeros(3*nb,5*nb),jacr';...
                                    zeros(4*nb,3*nb), bigJ, bigP',jacp';...
                                    zeros(nb,3*nb), bigP, zeros(nb,nb), zeros(nb,nc);...
                                    jacr, jacp, zeros(nc,nb+nc)];
                                obj = obj.assembleFullNRIJac();
                            end
                        case{'Full'}
                            obj.myPSY = [ bigM, zeros(3*nb,5*nb),jacr';...
                                zeros(4*nb,3*nb), bigJ, bigP',jacp';...
                                zeros(nb,3*nb), bigP, zeros(nb,nb), zeros(nb,nc);...
                                jacr, jacp, zeros(nc,nb+nc)];
                            obj = obj.assembleFullNRIJac();
                    end
                    %Step 2: Compute residual of system g(z)
                    
                    z = [rdd; pdd; lambP; lamb];
                    gtop = [ bigM, zeros(3*nb,5*nb),jac(1:end-nb,rindx)';...
                        zeros(4*nb,3*nb), bigJ, bigP',jac(1:end-nb,pindx)']*z;
                    obj = obj.getPhiJac(t);
                    %myPhi is [Phi;Phi_pnorm], change order here to be in
                    %same form as procedure;
                    PHI = (1/obj.myBetaH)*[obj.myPhi(obj.myNKC+1:end);...
                        obj.myPhi(1:obj.myNKC)];
                    res = [gtop-F; PHI];
                    %Step 3: Solve BigPHI*delz=-g(z) to get delz.
                    %Step 4: z(v+1) = z(v) + delz(v)
                    %Step 5: v = v+1
                    %Step 6: If converged, recompute position and velocity and
                    %store. t = t+1, restart step 0.
                    
                    % z = [ rdd; pdd; phi_p; phi]
                    delz = obj.myPSY\-res;
                    z = z+delz;
                    obj.myAcceleration(rindx,t) = z(1:3*nb);
                    obj.myAcceleration(pindx,t) = z(3*nb+1:7*nb);
                    obj.myLambda(:,t) = z(7*nb+1:end);
                    obj.myFreact(1:3*nb,t) = -phiR'*lamb;
                    obj.myFreact(3*nb+1:end,t) = -phiP'*lamb;

                    
                    
                    %% DEBUG ONLY
%                     [delzm, indxm] = max(abs(delz));
%                     disp(delzm);
%                     disp(indxm);
                    %%
                    
                    
                    %Check if converged
                    if abs(delz)<=tol
                        if t == 2 %If 2nd timestep, use BDF1 approx
                            r = obj.myPosition(rindx,1)+...
                                h*obj.myVelocity(rindx,1)+h^2*rdd;
                            p = obj.myPosition(pindx,1)+...
                                h*obj.myVelocity(pindx,1)+h^2*pdd;
                            rd = obj.myVelocity(rindx,1)+h*rdd;
                            pd = obj.myVelocity(pindx,1)+h*pdd;
                        else %If time>2, use BDF3 or 6
                            r = Cr+obj.myBetaH*z(1:3*nb);
                            p = Cp+obj.myBetaH*z(3*nb+1:7*nb);
                            rd = Cdr+sqrt(obj.myBetaH)*z(1:3*nb);
                            pd = Cdp+sqrt(obj.myBetaH)*z(3*nb+1:7*nb);
                        end
                        %Save all variables
                        obj.myPosition(rindx,t) = r;
                        obj.myPosition(pindx,t) = p;
                        obj.myVelocity(rindx,t) = rd;
                        obj.myVelocity(pindx,t) = pd;
                        obj.myITS(1,t) = it;
                        break
                    end
                end
                obj.myTimeStep = obj.myTimeStep+1;
            end
        end
        
        %%
        
        function obj = getFunctions(obj)
            %(NEEDS WORK) Numerically estimate the time
            %derivatives of the Kin Driving constraints of the system.
            for i = 1:obj.myNKC
                h = 10^-2;
                delh = obj.myStepSize*h;
                t = 0:delh:(obj.myFinalTime+obj.myStepSize);
                funt = str2func(obj.myKinCon{i}.fun);
                funt = funt(t);
                funDt = diff(funt)/delh;
                funDDt = diff(funDt)/delh;
                obj.myFunTimes.funt(i,:) = funt(1:h^(-1):numel(funt(:))-h^(-1));
                obj.myFunTimes.funDt(i,:) = funDt(1:h^(-1):numel(funt(:))-h^(-1));
                obj.myFunTimes.funDDt(i,:) = funDDt(1:h^(-1):numel(funt(:))-h^(-1)); %%%%%%%
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
            
            obj.myFlags = [1,1,0,0,0];
            
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
                    
                    [PHI,JAC,~,~,~] = fh(obj.myKinCon{i},obj,obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myPhi(i) = PHI;
                    obj.myJac(i,7*b1-6:7*b1) = JAC(1:7);
                elseif strcmp(b1,'ground')
                    obj.myqi = [0;0;0;1;0;0;0;...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [zeros(7,1);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [PHI,JAC,~,~,~] = fh(obj.myKinCon{i},obj,obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myPhi(i) = PHI;
                    obj.myJac(i,7*b2-6:7*b2) = JAC(8:14);
                else
                    obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [PHI,JAC,~,~,~] = fh(obj.myKinCon{i},obj,...
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
            obj.myFlags = [0,1,1,0,0];
            
            for i = 1:obj.myNKC
                b1 = obj.myKinCon{i}.body1;
                b2 = obj.myKinCon{i}.body2;
                fh = str2func(obj.myKinCon{i}.type);
                timeStep = obj.myTimeStep;
                %if body 2 is ground, ignore the outputs from the
                %constraint functions
                if strcmp(b2,'ground')
                    %                     obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                    %                         0;0;0;1;0;0;0];
                    if strcmp(obj.myKinCon{i}.type,'cons_cd')
                        obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                            obj.myKinCon{i}.sQjBAR' ;1;0;0;0];
                    else
                        obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                            0;0;0;1;0;0;0];
                    end
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        0;0;0;0;0;0;0];
                    [~,JAC,NU,~,~] = fh(obj.myKinCon{i},obj,obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myNu(i) = NU;
                    obj.myJac(i,7*b1-6:7*b1) = JAC(1:7);
                elseif strcmp(b1,'ground')
                    obj.myqi = [0;0;0;1;0;0;0;...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [zeros(7,1);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [~,JAC,NU,~,~] = fh(obj.myKinCon{i},obj,obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myNu(i) = NU;
                    obj.myJac(i,7*b2-6:7*b2) = JAC(8:14);
                    
                else
                    obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [~,JAC,NU,~,~] = fh(obj.myKinCon{i},obj,obj.myTimeStep,...
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
            obj.myFlags = [0,1,0,1,0];
            
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
                    %                     if strcmp(obj.myKinCon{i}.type,'cons_cd')
                    %                     obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                    %                         obj.myKinCon{i}.sQjBAR' ;1;0;0;0];
                    %                     else
                    %                     obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                    %                         0;0;0;1;0;0;0];
                    %                     end
                    
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        0;0;0;0;0;0;0];
                    [~,JAC,~,GAMMA,~] = fh(obj.myKinCon{i},obj,obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myGamma(i) = GAMMA;
                    obj.myJac(i,7*b1-6:7*b1) = JAC(1:7);
                elseif strcmp(b1,'ground')
                    obj.myqi = [0;0;0;1;0;0;0;...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [zeros(7,1);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [~,JAC,~,GAMMA,~] = fh(obj.myKinCon{i},obj,obj.myTimeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myGamma(i) = GAMMA;
                    obj.myJac(i,7*b2-6:7*b2) = JAC(8:14);
                    
                else
                    obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [~,JAC,~,GAMMA,~] = fh(obj.myKinCon{i},obj,obj.myTimeStep,...
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
                    % INDEX IS WRONG HERE FOR MULTIBODY??????? MB NOT
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
        function obj = getNRIjac(obj)
            obj.myFlags = [0,0,0,0,1];
            %             nb = obj.myNbodies;
            nc = obj.myNKC;
            timeStep = obj.myTimeStep;
            for i = 1:nc
                b1 = obj.myKinCon{i}.body1;
                b2 = obj.myKinCon{i}.body2;
                fh = str2func(obj.myKinCon{i}.type);
                if strcmp(b2,'ground')
                    obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        0;0;0;1;0;0;0];
                    
                    %                     if strcmp(obj.myKinCon{i}.type,'cons_cd')
                    %                         obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                    %                             obj.myKinCon{i}.sQjBAR' ;1;0;0;0];
                    %                     else
                    %                         obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                    %                             0;0;0;1;0;0;0];
                    %                     end
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        0;0;0;0;0;0;0];
                    [~,~,~,~,NRIjac] = fh(obj.myKinCon{i},obj,timeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    %If b2 is ground, the 0 entries of NRIjac will be added
                    %to the index of body1.
                    obj.myNRIjac{i} = struct('body1',b1,'body2',b1,'NRIjac',NRIjac);
                    
                elseif strcmp(b1,'ground')
                    obj.myqi = [0;0;0;1;0;0;0;...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [zeros(7,1);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [~,~,~,~,NRIjac] = fh(obj.myKinCon{i},obj,timeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myNRIjac{i} = struct('body1',b2,'body2',b2,'NRIjac',NRIjac);
                else
                    obj.myqi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    obj.myqdi = [obj.myVelocity(7*b1-6:7*b1,timeStep);...
                        obj.myVelocity(7*b2-6:7*b2,timeStep)];
                    [~,~,~,~,NRIjac] = fh(obj.myKinCon{i},obj,timeStep,...
                        obj.myFunTimes,obj.myqi,obj.myqdi,obj.myFlags);
                    obj.myNRIjac{i} = struct('body1',b1,'body2',b2,'NRIjac',NRIjac);
                    
                end
                
            end
        end
        %%
        function obj = assembleFullNRIJac(obj)
            nb = obj.myNbodies;
            nc = obj.myNKC;
            obj = obj.getNRIjac();
            t = obj.myTimeStep;
            %%%
            % CONTRIBUTION TO THE JACOBIAN FROM THE GENERALIZED FORCES NOT
            % IMPLEMENTED YET.
            %%%
            for i = 1:nc
                %Add the contribution of each constraint to the approriate
                %place in the global Jacobian matrix.
                %Slide 6 on:
                %http://sbel.wisc.edu/Courses/ME751/2016/Documents/lecture1019.pdf
                %PSY11
                
                idx = [obj.myNRIjac{i}.body1 obj.myNRIjac{i}.body2];
                %idx of the bodies involved in the constraint.
                for j =1:2 %row #
                    for k = 1:2 %Column #
                        %11
                        obj.myPSY(3*idx(j)-2:3*idx(j),3*idx(k)-2:3*idx(k)) = ...
                            obj.myPSY(3*idx(j)-2:3*idx(j),3*idx(k)-2:3*idx(k))+...
                            obj.myNRIjac{i}.NRIjac(3*j-2:3*j,3*k-2:3*k)*obj.myBetaH;
                        %21
                        obj.myPSY(3*nb+4*idx(j)-3:3*nb+4*idx(j),3*idx(k)-2:3*idx(k)) = ...
                            obj.myPSY(3*nb+4*idx(j)-3:3*nb+4*idx(j),3*idx(k)-2:3*idx(k)) +...
                            obj.myNRIjac{i}.NRIjac(6+4*j-3:6+4*j,3*k-2:3*k)*obj.myBetaH;
                        %12
                        obj.myPSY(3*idx(j)-2:3*idx(j),3*nb+4*idx(k)-3:3*nb+4*idx(k)) = ...
                            obj.myPSY(3*idx(j)-2:3*idx(j),3*nb+4*idx(k)-3:3*nb+4*idx(k)) +...
                            obj.myNRIjac{i}.NRIjac(3*j-2:3*j,6+4*k-3:6+4*k)*obj.myBetaH;
                        %22
                        obj.myPSY(3*nb+4*idx(j)-3:3*nb+4*idx(j),3*nb+4*idx(k)-3:3*nb+4*idx(k)) = ...
                            obj.myPSY(3*nb+4*idx(j)-3:3*nb+4*idx(j),3*nb+4*idx(k)-3:3*nb+4*idx(k)) +...
                            obj.myNRIjac{i}.NRIjac(6+4*j-3:6+4*j,6+4*k-3:6+4*k)*obj.myBetaH;
                        
                    end
                end
            end
            lam = obj.myLambda(1:nb,t);
            bigLam = diag(reshape(repmat(lam,1,4)',4*obj.myNbodies,1));
            for i = 1:nb
                p = obj.myPosition(7*i-3:7*i,t);
                pd = obj.myVelocity(7*i-3:7*i,t);
                pdd = obj.myAcceleration(7*i-3:7*i,t);
                J = diag(obj.myJ(3*i-2:3*i));
                G = obj.getG(p);
                Gd = obj.getG(pd);
                Gdd = obj.getG(pdd);
                a = J*G*pdd;
                aTIL = [ 0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
                T = [ 0 -a'; a -aTIL];
                touP = 8*Gd'*J*Gd;
                a2 = -J*Gd*p;
                a2TIL = [ 0 -a2(3) a2(2); a2(3) 0 -a2(1); -a2(2) a2(1) 0];
                T2 = [0 -a2'; a2 -a2TIL];
                touPd = 8*(T2-Gd'*J*G);
                BL = bigLam(4*i-3:4*i,4*i-3:4*i);
                obj.myPSY(3*nb+4*i-3:3*nb+4*i,3*nb+4*i-3:3*nb+4*i) = ...
                    obj.myPSY(3*nb+4*i-3:3*nb+4*i,3*nb+4*i-3:3*nb+4*i) +...
                    (4*(T-G'*J*Gdd)+BL-touP)*obj.myBetaH-...
                    sqrt(obj.myBetaH)*touPd;
                
                
            end
        end
        function obj = numericalJacobian(obj,timeStep)
            %Compute Jacobian of constraint equations via finite difference.
            flags = [ 1,0,0,0,0];
            nc = obj.myNKC;
            %             nb = obj.myNbodies;
            del = 10^-6;
            for i = 1:nc %cycle through all constraints
                b1 = obj.myKinCon{i}.body1;
                b2 = obj.myKinCon{i}.body2;
                fh = str2func(obj.myKinCon{i}.type);
                if strcmp(b2,'ground')
                    qi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        0;0;0;1;0;0;0];
                    n = 7;%See how many finite differences to compute
                else
                    qi = [obj.myPosition(7*b1-6:7*b1,timeStep);...
                        obj.myPosition(7*b2-6:7*b2,timeStep)];
                    n = 14;
                end
                k = 1; %used below to place jac items in correct spots
                
                for j = 1:n %cycle through each coordinate
                    qi2 = qi;
                    qi2(j) = qi2(j)+del; %perturb coordinate
                    [phi1,~,~,~,~] = fh(obj.myKinCon{i},obj,timeStep,...
                        obj.myFunTimes,qi,obj.myqdi,flags);
                    [phi2,~,~,~,~] = fh(obj.myKinCon{i},obj,timeStep,...
                        obj.myFunTimes,qi2,obj.myqdi,flags);
                    q = (phi1-phi2)/del;
                    
                    if j<=7
                        obj.myNumJac(i,7*b1-7+j) = q;
                    else
                        obj.myNumJac(i,7*b2-7+k) = q;
                        k = k+1;
                    end
                end
                
            end
        end
        function [ftsda] = getFTSDA(obj)
            %Haugh's book pg 446-447 (11.4)
            nb = obj.myNbodies;
            ftsda = zeros(7*nb,1);
            ts = obj.myTimeStep;
            for i = 1:obj.myNTSDA
                b1 = obj.myTSDA(i).body1;
                b2 = obj.myTSDA(i).body2;
                k = obj.myTSDA(i).k;
                l0 = obj.myTSDA(i).l0;
                c = obj.myTSDA(i).c;
%                 fun = obj.myTSDA(i).fun;
                
                
                b1idx = 7*b1-6:7*b1;
                ri = obj.myPosition(b1idx(1:3),ts);
                pi = obj.myPosition(b1idx(4:7),ts);
                rdi = obj.myVelocity(b1idx(1:3),ts);
                pdi = obj.myVelocity(b1idx(4:7),ts);
                
                if strcmp(b2,'ground')
                    rj = [0; 0; 0];
                    pj = [1; 0; 0; 0];
                    rdj = [0; 0; 0];
                    pdj = [0; 0; 0; 0];
                else
                    b2idx = 7*b2-6:7*b2;
                    rj = obj.myPosition(b2idx(1:3),ts);
                    pj = obj.myPosition(b2idx(4:7),ts);
                    rdj = obj.myVelocity(b2idx(1:3),ts);
                    pdj = obj.myVelocity(b2idx(4:7),ts);
                end
                Ai = obj.getA(pi);
                Aj = obj.getA(pj);
                dij = rj+Aj*obj.myTSDA(i).sQjBAR'-ri-Ai*obj.myTSDA(i).sPiBAR';
                l = norm(dij);
                
                sQjTIL = obj.getTIL(obj.myTSDA(i).sQjBAR);
                sPiTIL = obj.getTIL(obj.myTSDA(i).sPiBAR);
                Gj = obj.getG(pj);
                Gi = obj.getG(pi);
                
                ldot = (dij/l)'*(rdj-2*Aj*sQjTIL*Gj*pdj-rdi+2*Ai*sPiTIL*Gi*pdi);
                %ACTUATOR FUNCTION NOT IMPLEMENTED YET!!!!
                f = k*(l-l0)+c*ldot;
                
                if strcmp(b2,'ground')
                    ftsda(7*b1-6:7*b1) = ftsda(7*b1-6:7*b1)+(f/l)*[dij; 2*Gi'*sPiTIL*Ai'*dij];
                else
                    ftsda(7*b1-6:7*b1) = ftsda(7*b1-6:7*b1)+(f/l)*[dij; 2*Gi'*sPiTIL*Ai'*dij];
                    ftsda(7*b2-6:7*b2) = ftsda(7*b2-6:7*b2)-(f/l)*[dij; 2*Gj'*sQjTIL*Aj'*dij];
                    
                end
                
                
            end
        end
    end
    %found in Key Kinematic Equations:
    % http://sbel.wisc.edu/Courses/ME751/2010/Documents/kinematicsKeyFormulas.pdf
    methods(Static)
        function [G] = getG(p)
            e = p(2:end);
            eTIL = [ 0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];
            G = [-e, -eTIL+p(1)*eye(3)];
        end
        function [A] = getA(p)
            
            A = [ p(1).^2+p(2).^2-p(3).^2-p(4).^2 , 2*(p(2)*p(3)-p(1)*p(4)),...
                2*(p(2)*p(4)+p(1)*p(3)); 2*(p(2)*p(3)+p(1)*p(4)) ...
                p(1).^2-p(2).^2+p(3).^2-p(4).^2 , 2*(p(3)*p(4)-p(1)*p(2));...
                2*(p(2)*p(4)-p(1)*p(3)) , 2*(p(3)*p(4)+p(1)*p(2)), p(1).^2-p(2).^2-p(3).^2+p(4).^2];
        end
        function [til] = getTIL(a)
            til = [ 0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
        end
        
    end
end