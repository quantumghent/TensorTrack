classdef InfQP
    % Infinite Quasi-Particle states
    
    
    %% Properties
    properties
        mpsleft  UniformMps
        mpsright UniformMps
        X
        B
        VL
        p
    end
    
    properties (Dependent)
        AL
        AR
    end
    
    %% Constructors
    methods
        function qp = InfQP(varargin)
            if nargin == 0, return; end
            
            qp.mpsleft  = varargin{1};
            qp.mpsright = varargin{2};
            qp.X    = varargin{3};
            qp.VL   = varargin{4};
            qp.B    = varargin{5};
            qp.p    = varargin{6};
        end
    end
    
    methods (Static)
        function qp = new(fun, mpsleft, mpsright, p, charge)
            arguments
                fun                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                mpsleft
                mpsright = []
                p = 0
                charge = []
            end
            
            if isempty(mpsright), mpsright = mpsleft; end
            assert(period(mpsleft) == period(mpsright));
            
            dims = struct;
            dims.charges = charge;
            dims.degeneracies = ones(size(charge));
            
            AL = mpsleft.AL;
            for i = period(mpsleft):-1:1
                VL{i} = leftnull(AL{i});
                rVspace = rightvspace(VL{i});
                lVspace = leftvspace(mpsright, next(i, period(mpsleft)));
                if isempty(charge)
                    aspace = one(rVspace);
                else
                    aspace = rVspace.new(dims, false);
                end
                X{i} = Tensor.new(fun, rVspace', [aspace lVspace]);
            end
            
            qp = InfQP(mpsleft, mpsright, X, VL, [], p);
            qp.B = computeB(qp);
        end
        
        function qp = randnc(varargin)
            qp = InfQP.new(@randnc, varargin{:});
        end
    end
    
    
    %% Derived Properties
    methods
        function s = auxspace(qp, i)
            s = space(qp.X{i}, 3);
        end
        
        function al = get.AL(qp)
            al = qp.mpsleft.AL;
        end
        
        function ar = get.AR(qp)
            ar = qp.mpsright.AR;
        end
        
        function B = computeB(qp)
            for w = period(qp):-1:1
                B{w} = multiplyright(qp.VL{w}, qp.X{w});
            end
        end
        
        function X = computeX(qp)
            for w = period(qp):-1:1
                X{w} = tracetransfer(qp.VL{w}', qp.B{w});
            end
        end
        
        function bool = istrivial(qp)
            bool = qp.p == 0 && istrivial(auxspace(qp, 1));
        end
        
        function rho = fixedpoint(qp, type, w)
            % compute the fixed point of the transfer matrix of an mps.
            %
            % Usage
            % -----
            % :code:`rho = fixedpoint(mps, type, w)`
            %
            % Arguments
            % ---------
            % mps : :class:`UniformMps`
            %   input state.
            %
            % type : char
            %   specification of the type of transfer matrix:
            %   general format: sprintf(%c_%c%c, side, top, bot) where side is 'l' or 'r' to
            %   determine which fixedpoint, and top and bot are 'L' or 'R' to specify
            %   whether to use AL or AR in the transfer matrix.
            %
            % w : integer
            %   position within the mps unitcell of the fixed point.
            
            arguments
                qp
                type {mustBeMember(type, ...
                    {'l_LL' 'l_LR' 'l_RL' 'l_RR' 'r_LL' 'r_LR' 'r_RL' 'r_RR'})}
                w = strcmp(type(1), 'l') * 1 + strcmp(type(1), 'r') * period(qp)
            end
            
            switch type(1)
                case 'l'
                    rho = fixedpoint(qp.mpsleft, type, w);
                case 'r'
                    rho = fixedpoint(qp.mpsright, type, w);
                otherwise
                    error('invalid type');
            end
        end
    end
    
    methods
        function p = period(qp)
            p = length(qp.X);
        end
        
        function type = underlyingType(qp)
            type = underlyingType(qp.X{1});
        end
    end
end

