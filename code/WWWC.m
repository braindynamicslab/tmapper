classdef WWWC
%WWWC The Wong-Wang-Wilson-Cowan hybrid model, which takes the form of
% ----------------------------------------------------------------------------------------
%   dS_i,E/dt = - S_i,E/tau_E + (1-S_i,E) * gamma_E * H_E(w_EE.*S_i,E - w_IE.*S_i,I + I_E + I_i,G)
%   dS_i,I/dt = - S_i,I/tau_I + (1-S_i,I) * gamma_I * H_I(w_EI.*S_i,E - w_II.*S_i,I + I_I)
% where the global input
%   I_i,G = G * sum_j [C_ij*S_i,E] where C_ij = 0 for i=j
% and the transfer function
%   H_p(x) = (r_max + (a_p*x-b_p-r_max) / (1-exp(d_p * (a_p*x - b_p -
%   r_max))) )/ (1-exp(-d_p*(a_p*x-b_p)))
% for p=E or I.
% ----------------------------------------------------------------------------------------
% The model can be thought of as a high-dimensional generalization of
% Wong-Wang model which is more consistent with the form of the
% Wilson-Cowan model. Importantly, the model use a sigmoidal transfer
% function that matches the linear threshold transfer function as used in
% Wong-Wang (2006)/Deco et al (2014JON) asymptotically for lower level of
% of activity. Note that sigmoidal transfer function will cap the firing
% rate such that the model is more realistic. 
% 
%   obj = WWWC(C)   create a model based on a connectivity matrix C
%   obj = WWWC(C, G)    create a model based on a connectivity matrix C and
%                       a global coupling parameter G (scalar, or scalar
%                       function of time)
%   obj.EularSolver(X_0)    solve system numerically using Eular-Maruyama
%                           method. 
%   obj.HeunSolver(X_0)     solve system numerically using Heun method. 
%{
created by MZ, 7-2-2019
%}
    
    properties
        N; % number of nodes 
        C; % connectivity matrix with zero diagonal
        
        % --- other model parameters
        tau_E = 0.1; % time constant for NMDA
        tau_I = 0.01; % time constant for GABA
        gamma_E = 0.641; % kinetic parameter for excitatory popu
        gamma_I = 1; % kinetic parameter for inhibitory popu
        w_EE = 0.25; % excitatory-to-excitatory coupling (scalar or N-element vector)
        w_IE = 0.25; % inhibitory-to-excitatory coupling (scalar or N-element vector)
        w_EI = 0.4; % excitatory-to-inhibitory coupling (scalar or N-element vector)
        w_II = 0.05; % inhibitory-to-inhibitory coupling (scalar or N-element vector)
        I_E = 0; % ambient input to excitatory popu
        I_I = 0.1; % ambient input to inhibitory popu
        G = @(t) 1; % global coupling (can be a function of time)
        % ~ parameters of the transfer function
        a_E = 310; 
        b_E = 125;
        d_E = 0.16;
        a_I = 615;
        b_I = 177;
        d_I = 0.087;
        r_max = 500; % maximal firing rate
        % ~ noise term
        sigma = 0;
        
        % --- simulation parameters
        X_0;
        dt = 0.01; %(s) step size
        T = 1; % %(s) total integration time
        t; % (s) time axis
        Nt; % (count) number of time samples
        
        % simulation results
        X; 
    end
    
    methods
        function obj = WWWC(C, G)
            %WWWC Construct an instance of this class
            %   
            if nargin>1
                obj.G = G;
            end
            obj.C = C;        
        end
        
        %% ===== the dynamical system ===== %%
        function dX = D(obj,t,X)
            %D the differential operator
            % input:
            %   X: can be vector of 2N elements. The first N elements are states
            %   of the excitatory population; the last N elements are
            %   states of the inhibitory population.
            dX=nan(size(X));
            SE=X(1:obj.N,:);
            SI=X(obj.N+1:end,:);
            
            dX(1:obj.N,:) = - SE/obj.tau_E + (1-SE)*obj.gamma_E.*...
                obj.H_E(obj.w_EE.*SE - obj.w_IE.*SI + obj.I_E + obj.G(t)*obj.C*SE);
            dX(obj.N+1:end,:) = -SI/obj.tau_I + (1-SI)*obj.gamma_I.*...
                obj.H_I(obj.w_EI.*SE - obj.w_II.*SI + obj.I_I);
        end
        
        function r = H_E(obj,x)
            %H_E the excitatory transfer function
            r = (obj.r_max+(obj.a_E*x-obj.b_E-obj.r_max)./(1-exp(obj.d_E*(obj.a_E*x-obj.b_E-obj.r_max))))...
                ./(1-exp(-obj.d_E*(obj.a_E*x-obj.b_E)));
        end
        
        function r = H_I(obj,x)
            %H_I the inhibitory transfer function
            r = (obj.r_max+(obj.a_I*x-obj.b_I-obj.r_max)./(1-exp(obj.d_I*(obj.a_I*x-obj.b_I-obj.r_max))))...
                ./(1-exp(-obj.d_I*(obj.a_I*x-obj.b_I)));
        end
        
        %% ===== methods of simulations ===== %%
        function obj = EularSolver(obj,X_0)
            %EULARSOLVER solve equation with Eular-Maruyama method
            if nargin>1 && ~isempty(X_0)
                obj.X_0 = X_0;
            end
            
            % initialize
            obj.t=(0:obj.dt:obj.T)';
            obj.Nt=length(obj.t);
            sln=nan(2*obj.N,obj.Nt);
            sln(:,1)=obj.X_0;
            
            % simulate
            for tidx=2:obj.Nt
                sln(:,tidx) = sln(:,tidx-1) ...
                            + obj.dt*obj.D(obj.t(tidx),sln(:,tidx-1)) ... drift term
                            + obj.sigma*sqrt(obj.dt)*randn(2*obj.N,1); % diffusion term
            end
            obj.X = sln';
        end
        
        function obj = HeunSolver(obj,X_0)
            % solve euqation with Heun Stochastic Scheme (Mannella 2002)
            if nargin>1 && ~isempty(X_0)
                obj.X_0 = X_0;
            end
            
            % initialize
            obj.t=(0:obj.dt:obj.T)';
            obj.Nt=length(obj.t);
            sln=nan(2*obj.N,obj.Nt);
            sln(:,1)=obj.X_0;
            
            % simulate
            for tidx=2:obj.Nt
                noise = randn(2*obj.N,1);
                dX = obj.D(obj.t(tidx),sln(:,tidx-1));
                % --- step 1: find x_{i+1} using Eular
                slntmp = sln(:,tidx-1) ...
                            + obj.dt * dX ... drift term (first order approx)
                            + obj.sigma*sqrt(obj.dt)*noise; % diffusion term
                
                % --- step 2: average now and future derivative
                sln(:,tidx) = sln(:,tidx-1) ...
                            + obj.dt/2*(dX + obj.D(obj.t(tidx),slntmp)) ... drift term (second order approx)
                            + obj.sigma*sqrt(obj.dt)*noise; % diffusion term
            end
            obj.X = sln';
        end
        
        %% ===== setting parameters ===== %%
        function obj = set.C(obj, C)
            obj.C = zerodiag(C);  
            obj.N = size(obj.C,1);
        end
        
        function obj = set.G(obj, G)
            if isa(G,'function_handle')
                obj.G = G;
            elseif isnumeric(G)
                obj.G = @(t) G;
            else
                error('G must be a function handle or a number.')
            end
        end
    end
end

