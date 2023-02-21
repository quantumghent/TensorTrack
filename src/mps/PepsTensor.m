classdef PepsTensor
    % Generic PEPS tensor objects that have a notion of virtual and physical legs.
    
    properties
        var
    end
    
    
    %% Constructors
    methods
        function A = PepsTensor(tensor)
            arguments
                tensor = []
            end
            
            if iscell(tensor)
                for i = length(tensor):-1:1
                    A(i) = PepsTensor(tensor{i});
                end
                return
            end
            
            if ~isempty(tensor)
                A.var = tensor;
            end
        end
    end
    
    
    %% Properties
    methods
        function s = space(A, varargin)
            s = space(A.var, varargin{:});
        end
        
        function n = nspaces(A)
            n = nspaces(A.var);
        end
        
        function s = pspace(A)
            s = space(A, 1);
        end
        
        function s = westvspace(A)
            s = space(A, 2);
        end
        
        function s = southvspace(A)
            s = space(A, 3);
        end
        
        function s = eastvspace(A)
            s = space(A, 4);
        end

        function s = northvspace(A)
            s = space(A, 5);
        end
        
        function s = vspace(A)
            s = space(A, 2:5);
        end
        
        function cod = codomain(A)
            cod = A.var.codomain;
        end
        
        function dom = domain(A)
            dom = A.var.domain;
        end
        
        function r = rank(A)
            r = rank(A.var);
        end
    end
    
    
    %% Linear Algebra
    methods
        function d = dot(A, B)
            % TODO?
        end
        
        function A = repartition(A, varargin)
            A.var = repartition(A.var, varargin{:});
        end
        
        function A = tpermute(A, varargin)
            A.var = tpermute(A.var, varargin{:});
        end
        
        function A = plus(varargin)
            for i = 1:2
                if isa(varargin{i}, 'PepsTensor')
                    varargin{i} = varargin{i}.var;
                end
            end
            A = plus(varargin{:});
        end
        
        function A = minus(varargin)
            for i = 1:2
                if isa(varargin{i}, 'PepsTensor')
                    varargin{i} = varargin{i}.var;
                end
            end
            A = minus(varargin{:});
        end
        
        function n = norm(A)
            n = norm(A.var);
        end
        
        function A = conj(A)
            A.var = conj(A.var);
        end

        function A = twist(A, varargin)
            A.var = twist(A.var, varargin{:});
        end
                
        function t = ctranspose(t)
            % Compute the adjoint of a tensor. This is defined as swapping the codomain and
            % domain, while computing the adjoint of the matrix blocks.
            %
            % Usage
            % -----
            % :code:`t = ctranspose(t)`
            % :code:`t = t'`
            %
            % Arguments
            % ---------
            % t : :class:`Tensor`
            %   input tensor.
            %
            % Returns
            % -------
            % t : :class:`Tensor`
            %   adjoint tensor.
            
            t.var = t.var';
            
            t = permute(t, ndims(t):-1:1);
        end
        
        function C = tensorprod(varargin)
            for i = 1:2
                if isa(varargin{i}, 'PepsTensor')
                    varargin{i} = varargin{i}.var;
                end
            end
            C = tensorprod(varargin{:});
        end
        
        function type = underlyingType(A)
            type = underlyingType(A.var);
        end
    end
    
    
    %% Converters
    methods
        function t = Tensor(A)
            t = full(A.var);
        end
        
        function t = SparseTensor(A)
            t = reshape([A.var], size(A));
            t = sparse(t);
        end
    end
end

