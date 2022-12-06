classdef BallWindBOLD
    %BALLWINDBOLD Balloon-Windkessel model of BOLD dynamics
    %   Following Friston et al (2003).
    % obj = BallWindBOLD(Z)
    %   convert neurodynamic function Z to BOLD activity. Here Z is
    %   Nt-by-Nnodes matrix, where Nt is the number of time samples, and
    %   Nnodes is the number of brain areas involved. The BOLD dynamics is
    %   returned as a property obj.B. 
    %   The model of BOLD for a single brain area takes the form:
    %   -----------------------------------------
    %   ds/dt = z - kappa*s - gamma * (f-1)
    %   df/dt = s
    %   tau*dv/dt = f - v^(1/alpha)
    %   tau*dq/dt = f/rho*[1- (1-rho)^(1/f)] - v^(1/alpha - 1)*q
    %   B = V0*[k1*(1-q) + k2*(1-q/v)+k3*(1-v)]
    %   -----------------------------------------
    %   Here [s,f,v,q] are internal state variables generating the BOLD
    %   signal. B is the resulted BOLD signal. Here the initial condition
    %   of the internal state variables are [0, 1, 1, 1], which is at
    %   equilibrium without external input. See reference below for more
    %   details:
    %       Friston, K. J., Harrison, L., & Penny, W. (2003). Dynamic
    %       causal modelling. NeuroImage, 19(4), 1273?1302.
    %       http://doi.org/10.1016/S1053-8119(03)00202-7
    % obj = obj.EularSolver
    %   simulate BOLD signal. Here two parameters (dt & T) need to match
    %   that of the neural dynamic input Z. The class automatically impose
    %   a constraint on the relation between dt and T, i.e. dt*(Nt-1)=T.
    %   Here Nt is fixed from input data, should not be altered.
%{
created by MZ, 4/9/2018 
(mod 1) redefined obj.t everytime there's a change in dt or T (7-12-2019)
%}
    
    
    properties
        kappa = 0.65; % (s^-1) rate of signal decay
        gamma = 0.41; % (s^-1) rate of flow-dependent elimination
        tau = 0.98; %(s) hemodynamic transit time
        alpha = 0.32; % Grubb's exponent
        rho = 0.34; % resting oxygen extraction fraction 
        V0 = 0.02;
        k1 = 7*0.34; % k1=7*rho
        k2 = 2;
        k3 = 2*0.34-0.2; % k3 = 2*rho -0.2
        % --- numeric paramters
        dt = 0.001;
        T = 1;
        X0 = [0, 1, 1, 1];
        Nnodes; % number of brain areas
        Nvar = 4; % number of dynamic variables
        Nt; % number of samples
        t;
        
        % --- solutions
        Z; % brain activity Nt-by-Nnodes
        X; % dynamic variables
        B; % bold signal
        
        % --- calculate connectivity
        t_transient = 0; % duration of transient, to be removed from calulation of functional connectivity
        FC; % functional connectivity
        FCz; % Fisher-z transformation of FC = atanh(FC) -- ** a vector of subdiagonal elements **
        
    end
    
    methods
        function obj = BallWindBOLD(Z)
            %BALLWINDBOLD Construct an instance of this class
            
            obj.Z = Z;
            [obj.Nt,obj.Nnodes]=size(Z);
            obj.X = nan(obj.Nnodes,obj.Nvar,obj.Nt); % dynamic states
            obj.B = nan(obj.Nt, obj.Nnodes); % bold signal
            obj.T = (obj.Nt-1)*obj.dt;
            obj.t = (0:obj.Nt-1)'*obj.dt;
        end
        
        function dX = Fun(obj,nt,X)
            %Fun dynamic equation takes neuronal activity z as input
            %   input:
%                   X: state variables [z, s, f, v, q] X Nnodes rows.
            dX = nan(size(X));
            % c1:s
            dX(:,1) = obj.Z(nt,:)' - obj.kappa*X(:,1)-obj.gamma*(X(:,2)-1);
            % c2:f
            dX(:,2) = X(:,1);
            % c3:v
            dX(:,3) =(X(:,2)-X(:,3).^(1/obj.alpha))/obj.tau;
            % c4:q
            dX(:,4) = (X(:,2).*(1-(1-obj.rho).^(1./X(:,2)))/obj.rho - X(:,3).^(1/obj.alpha-1).*X(:,4))/obj.tau;
            
        end
        
        function obj = EularSolver(obj)
            % hemodynamic state variables
            obj.X(:,:,1) = repmat(obj.X0,obj.Nnodes,1);
            for n = 2:obj.Nt
                obj.X(:,:,n) = obj.X(:,:,n-1) + obj.dt*obj.Fun(n-1,obj.X(:,:,n-1));
            end
            % BOLD
            obj.B = squeeze(obj.V0*(obj.k1*(1-obj.X(:,4,:)) + obj.k2*(1-obj.X(:,4,:)./obj.X(:,3,:)) + obj.k3*(1-obj.X(:,3,:))))';
        end
        
        function FC = get.FC(obj)
            % functional connectivity
            FC = corrcoef(obj.B(obj.t>obj.t_transient,:));
        end
        
        function FCz = get.FCz(obj)
            FCz = fisherzFC(obj.FC);
        end
        
        function obj = set.dt(obj, dt)
            obj.dt = dt;
            if obj.T ~= dt*(obj.Nt-1)
                obj.T = dt*(obj.Nt-1);
                obj.t = (0:obj.Nt-1)'*dt;
            end
            disp(['set dt=' num2str(obj.dt) ' and T=' num2str(obj.T)])
        end
        
        function obj = set.T(obj, T)
            obj.T = T;
            if obj.dt ~= T/(obj.Nt-1)
                obj.dt = T/(obj.Nt-1);
                obj.t = (0:obj.Nt-1)'*obj.dt;
            end
            disp(['set dt=' num2str(obj.dt) ' and T=' num2str(obj.T)])
        end
        
        function obj = set.rho(obj,rho)
            % set "rho", which will change the value of k1 and k3
            obj.rho = rho;
            obj.k1 = 7*rho;
            obj.k3 = 2*rho-0.2;
            disp(['set rho=' num2str(obj.rho) ', k1=' num2str(obj.k1), ', k2=' num2str(obj.k2), ', k3=' num2str(obj.k3)])
            
        end
    end
end

