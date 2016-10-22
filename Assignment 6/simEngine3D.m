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
        %Number of constraint equations (7*nbodies), 6DOF + 1 Euler Parameter normalization.
        myNACE %Number of algebraic constraint equations
        %(=myNKC if kinematic analysis)
        
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
        myJac %
        myPhi %
        myNu %
        myGamma %
        myFlags % 4x1
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
            obj.myNKC = numel(obj.myProblemID.constraints);
            obj.myNACE = 7 * obj.myNbodies;
            %Initialize position velocity and accel of bodies to 0.
            obj.myPosition = zeros(7*obj.myNbodies,...
                obj.myProblemInfo.outputSteps);
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
            
        end
        function obj = getFunctions(obj)
            %Really bad implementation to numerically estimate the time
            %derivatives of the Kin Driving constraints of the system.
            for i = 1:obj.myNKC
                h = 10^-4;
                delh = obj.myStepSize*h;
                t = 0:delh:obj.myFinalTime;
                funt = str2func(obj.myKinCon{i}.fun);
                funt = funt(t);
                funDt = diff(funt)/delh;
                funDDt = diff(funDt)/delh;
                obj.myFunTimes.funt(i,:) = funt(1:h^(-1):numel(funt(:)));
                obj.myFunTimes.funDt(i,2:end) = funDt(1:h^(-1):numel(funt(:))-1);
                obj.myFunTimes.funDDt(i,2:end) = funDDt(1:h^(-1):numel(funt(:))-2);
                
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
                %                 disp(obj.myCurrentTime); %Uncomment to display time of
                %                 iteration
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
                obj.myPhi(6*obj.myNbodies+i) = p'*p-1;
                obj.myJac(6*obj.myNbodies+i,7*i-3:7*i) = 2*p';
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
            obj.myGamma = zeros(obj.myNACE,1);
            obj.myTimeStep = 1;
            obj.myCurrentTime = 0 + obj.myStepSize;
            obj.myFlags = [0,1,0,1];
            
            for timestep = 1:obj.myFinalTime/obj.myStepSize
                obj = getGammaJac(obj,timestep);
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
                timeStep = obj.myTimeStep;
                
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
                obj.myGamma(6*obj.myNbodies+i) = 2*(pdi'*pdi);
                obj.myJac(6*obj.myNbodies+i,7*i-3:7*i) = 2*p';
            end
        end
        %%
    end
    
end