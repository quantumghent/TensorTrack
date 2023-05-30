classdef PepsTensor
    % Generic PEPS tensor object that hos a notion of virtual and physical legs.
    %   This object represents the PEPS tensor at a single site as a rank (1, 4) tensor,
    %   where the physical index lies in the codomain and the virtual indices lie in the
    %   domain.
    %
    %                               5   1
    %                               |  /
    %                               v ^
    %                               |/
    %                        2 ->-- O --<- 4
    %                               |
    %                               ^
    %                               |
    %                               3
    properties
        var
    end
    
    
    %% Constructors
    methods
        function A = PepsTensor(tensor)
            arguments
                tensor = []
            end
                        
            if ~isempty(tensor)
                A.var = tensor;
            end
        end
    end
    
    
    %% Properties
    methods
        function n = nspaces(A)
            n = 5;
        end
        
        function s = space(A, varargin)
            s = space(A.var, varargin{:});
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
        function A = repartition(A, varargin)
            A.var = repartition(A.var, varargin{:});
        end
        
        function A = tpermute(A, varargin)
            A.var = tpermute(A.var, varargin{:});
        end
        
        function A = plus(A, B)
            arguments
                A PepsTensor
                B PepsTensor
            end
            A.var = plus(A.var, B.var);
        end
        
        function A = minus(A, B)
            arguments
                A PepsTensor
                B PepsTensor
            end
            A.var = minus(A.var, B.var);
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
            t.var = t.var';
            t = permute(t, ndims(t):-1:1);
        end

        function t = dagger(t)
            leftflip = Tensor.eye(westvspace(t), eastvspace(t));
            rightflip = twist(leftflip', 2);
            pflip =  Tensor.eye(pspace(t), pspace(t)');
            t.var = tpermute(conj(t.var), [1,2,5,4,3]);
            t.var = contract(t.var, [1,2,-3,3,-5], pflip, [1,-1], leftflip, [-2,2], rightflip, [3,-4]);
        end

        function t = rot90(t)
            t = tpermute(t, [1, 3, 4, 5, 2], rank(t));
        end
        
        function t = rot270(t)
            t = tpermute(t, [1, 5, 2, 3, 4], rank(t));
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
    

    %% Static constructors
    methods (Static)
        function t = new(fun, pspace, westvspace, southvspace, eastvspace, northvspace)
            arguments
                fun
                pspace
                westvspace
                southvspace
                eastvspace = westvspace'
                northvspace = southvspace'
            end
            
            t = PepsTensor(Tensor.new(fun, pspace, [westvspace, southvspace, eastvspace, northvspace]'));
        end
        
        function t = rand(varargin)
            t = PepsTensor.new(@rand, varargin{:});
        end
        
        function t = randn(varargin)
            t = PepsTensor.new(@randn, varargin{:});
        end
        
        function t = randc(varargin)
            t = PepsTensor.new(@randc, varargin{:});
        end
        
        function t = randnc(varargin)
            t = PepsTensor.new(@randnc, varargin{:});
        end
    end
end

