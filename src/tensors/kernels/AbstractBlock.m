classdef (Abstract) AbstractBlock
    % Abstract structure for storing tensor data.
    %
    %   This represents the blocks in the block-diagonal decomposition of a general tensor.
    
    %#ok<*INUSD>
    %#ok<*INUSL>
    %#ok<*STOUT>
    
    %% Constructor
    methods (Static)
        function X = new(codomain, domain)
            % Dispatch block type based on codomain and domain.
            %
            % Arguments
            % ---------
            % fun : :class:`function_handle`
            %   initialising function for the tensor data, with signature
            %   :code:`data = fun(dims)` where dims is a row vector of dimensions.
            %
            % codomain : :class:`.AbstractSpace`
            %   vector of vector spaces that form the codomain.
            %
            % domain : :class:`.AbstractSpace`
            %   vector of vector spaces that form the domain.
            %
            % Returns
            % -------
            % X : :class:`.AbstractBlock`
            %   tensor data.
            
            if isa(codomain, 'CartesianSpace') || isa(codomain, 'ComplexSpace') || ...
                    isa(domain, 'CartesianSpace') || isa(domain, 'ComplexSpace')
                X = TrivialBlock(codomain, domain);
                return
            end
            
            assert(~isa(codomain, 'SumSpace') && ~isa(domain, 'SumSpace'), ...
                'tensors:argerror', 'Cannot construct tensor with SumSpaces.');
            
            global cache
            if isempty(cache), cache = LRU; end
            
            if Options.CacheEnabled()
                key = GetMD5({codomain, domain}, 'Array', 'hex');
                med = get(cache, key);
                if isempty(med)
                    if braidingstyle(codomain, domain) == BraidingStyle.Abelian && ...
                            fusionstyle(codomain, domain) == FusionStyle.Unique
                        med = AbelianBlock(codomain, domain);
                    else
                        med = MatrixBlock(codomain, domain);
                    end
                    cache = set(cache, key, med);
                end
            else
                if braidingstyle(codomain, domain) == BraidingStyle.Abelian && ...
                        fusionstyle(codomain, domain) == FusionStyle.Unique
                    med = AbelianBlock(codomain, domain);
                else
                    med = MatrixBlock(codomain, domain);
                end
            end
            X = med;
        end
    end
    
    
    %% Required methods
    methods
        function Y = axpby(a, X, b, Y, p, map)
            % Compute :code:`Y = permute(X, p) .* a + Y .* b`.
            % This method is the computationally critical method of this class, thus has
            % special cases for scalar multiplication (:code:`a == 0`), addition
            % (:code:`nargin == 4`), and various optimizations when :code:`a == 1`,
            % :code:`b == 0 || b == 1`. Additionally, this method should not be called
            % directly as it should not perform any error checks.
            %
            % Arguments
            % ---------
            % a : :class:`double`
            %   scalar to multiply with X.
            %
            % X : :class:`.AbstractBlock`
            %   list of source blocks.
            %
            % b : :class:`double`
            %   scalar to multiply with Y.
            %
            % Y : :class:`.AbstractBlock`
            %   list of destination blocks.
            %
            % p : :class:`int`
            %   permutation vector for X.
            %
            % map : (sparse) :class:`double`
            %   coefficient matrix for permuting X.
            %
            % Returns
            % -------
            % Y : :class:`.AbstractBlock`
            %   Result of computing :code:`Y = permute(X, p) .* a + Y .* b`.
            %
            % Note
            % ----
            % This is an abstract method that should be overloaded for each subtype.
            error('This method should be overloaded.');
        end
        
        function [mblocks, mcharges] = matrixblocks(b)
            % Extract a list of coupled matrix blocks.
            %
            % Arguments
            % ---------
            % b : :class:`.AbstractBlock`
            %   list of input data.
            %
            % Returns
            % -------
            % mblocks : :class:`cell`
            %   list of non-zero coupled matrix blocks, sorted according to its charge.
            %
            % mcharges : :class:`.AbstractCharge`
            %   list of coupled charges.
            %
            % Note
            % ----
            % This is an abstract method that should be overloaded for each subtype.
            
            error('This method should be overloaded.');
        end
        
        function C = mul(C, A, B, a, b)
            % Compute :code:`C = (A .* a) * (B .* b)`.
            % Compute the matrix product of two source tensor structures, and store the
            % result in the destination tensor structure. This method should not perform any
            % error checks.
            %
            % Arguments
            % ---------
            % C : :class:`.AbstractBlock`
            %   location to store the result
            %
            % A : :class:`.AbstractBlock`
            %   first matrix factor
            %
            % B : :class:`.AbstractBlock`
            %   second matrix factor
            %
            % a : :class:`double` = 1
            %   first scalar factor
            %
            % b : :class:`double` = 1
            %   second scalar factor
            %
            % Returns
            % -------
            % C : :class:`.AbstractBlock`
            %   Result of computing C = (A .* a) * (B .* b)
            %
            % Note
            % ----
            % This is an abstract method that should be overloaded for each subtype.
            
            error('This method should be overloaded.');
        end
        
        function tblocks = tensorblocks(b)
            % Extract a list of uncoupled tensor blocks.
            %
            % Arguments
            % ---------
            % b : :class:`.AbstractBlock`
            %   list of input data.
            %
            % Returns
            % -------
            % tblocks : :class:`cell`
            %   list of non-zero uncoupled tensor blocks, sorted according to the coupled
            %   charge and then in column-major order according to the uncoupled charges.
            %
            % Note
            % ----
            % This is an abstract method that should be overloaded for each subtype.
            
            error('This method should be overloaded.');
        end
    end
    
    
    %% Optional methods
    methods
        function Y = axpy(a, X, Y, p, map)
            % Compute :code:`Y = permute(X, p) .* a + Y`.
            % This method is a convenience method that automatically falls back on
            % :code:`Y = axpby(a, X, 1, Y, p, map)`, but can be overloaded if the additional
            % efficiency is desired.
            %
            % Arguments
            % ---------
            % a : :class:`double`
            %   scalar to multiply with X.
            %
            % X : :class:`.AbstractBlock`
            %   list of source blocks.
            %
            % Y : :class:`.AbstractBlock`
            %   list of destination blocks.
            %
            % p : :class:`int`
            %   permutation vector for X.
            %
            % map : (sparse) :class:`double`
            %   coefficient matrix for permuting X.
            %
            % Returns
            % -------
            % Y : :class:`.AbstractBlock`
            %   Result of computing :code:`Y = permute(X, p) .* a + Y`.
            
            Y = axpby(a, X, 1, Y, p, map);
        end
        
        function Y = minus(X, Y)
            % Subtraction of :code:`X` and :code`Y`.
            %
            % Usage
            % -----
            % :code:`Y = minus(X, Y)`
            % 
            % :code:`Y = X - Y`
            %
            % Arguments
            % ---------
            % X : :class:`.AbstractBlock`
            %   first list of input data.
            %
            % Y : :class:`.AbstractBlock`
            %   second list of input data.
            %
            % Returns
            % -------
            % Y : :class:`.AbstractBlock`
            %   list of output data.
            
            Y = axpby(1, X, -1, Y);
        end
        
        function Y = plus(X, Y)
            % Addition of :code:`X` and :code:`Y`.
            %
            % Usage
            % -----
            % :code:`Y = plus(X, Y)`
            % 
            % :code:`Y = X + Y`
            %
            % Arguments
            % ---------
            % X : :class:`.AbstractBlock`
            %   first list of input data.
            %
            % Y : :class:`.AbstractBlock`
            %   second list of input data.
            %
            % Returns
            % -------
            % Y : :class:`.AbstractBlock`
            %   list of output data.
            
            Y = axpby(1, X, 1, Y);
        end
        
        function Y = times(Y, a)
            % Scalar multiplication of :code:`Y` and :code:`a`.
            %
            % Usage
            % -----
            % :code:`Y = times(Y, a)`
            % 
            % :code:`Y = Y .* a`
            %
            % Arguments
            % ---------
            % Y : :class:`.AbstractBlock`
            %   list of input data.
            %
            % a : :class:`double`
            %   scalar factor.
            %
            % Returns
            % -------
            % Y : :class:`.AbstractBlock`
            %   list of output data.
            
            Y = axpby(0, [], a, Y);
        end
        
        function Y = rdivide(Y, a)
            % Scalar division of Y and a.
            %
            % Usage
            % -----
            % :code:`Y = rdivide(Y, a)`
            % :code:`Y = Y ./ a`
            %
            % Arguments
            % ---------
            % Y : :class:`.AbstractBlock`
            %   list of input data.
            %
            % a : :class:`double`
            %   scalar factor.
            %
            % Returns
            % -------
            % Y : :class:`.AbstractBlock`
            %   list of output data.
            
            Y = axpby(0, [], 1/a, Y);
        end
        
        function A = uplus(A)
            % Unary plus. Equivalent to making a copy.
            %
            % Usage
            % -----
            % :code:`A = uplus(A)`
            % 
            % :code:`A = +A`
            %
            % Arguments
            % ---------
            % A : :class:`.AbstractBlock`
            %   list of input data.
            %
            % Returns
            % -------
            % A : :class:`.AbstractBlock`
            %   list of output data.
            
        end
        
        function Y = uminus(Y)
            % Unary minus. Computes the additive inverse.
            %
            % Usage
            % -----
            % :code:`A = uminus(A)`
            % 
            % :code:`A = -A`
            %
            % Arguments
            % ---------
            % A : :class:`.AbstractBlock`
            %   list of input data.
            %
            % Returns
            % -------
            % A : :class:`.AbstractBlock`
            %   list of output data.
            
            Y = Y .* (-1);
        end
        
        function t = underlyingType(X)
            % The scalar type of the tensor.
            %
            % Arguments
            % ---------
            % X : :class:`.AbstractBlock`
            %   input block.
            % 
            % Returns
            % -------
            % type : :class:`char`
            %   the scalar type of the data, which is one of the following:
            %
            %   - :class:`single` or :class:`double`
            %   - :class:`logical`
            %   - :class:`int8`, :class:`int16`, :class:`int32` or :class:`int64`
            %   - :class:`uint8`, :class:`uint16`, :class:`uint32` or :class:`uint64`
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end
        
        function X = ctranspose(X)
            % Adjoint of a tensor.
            %
            % Usage
            % -----
            % :code:`X = ctranspose(X)`
            %
            % :code:`X = X'`
            %
            % Arguments
            % ---------
            % X : :class:`.AbstractBlock`
            %   list of input data.
            %
            % Returns
            % -------
            % X : :class:`.AbstractBlock`
            %   list of adjoint output data.
            
        end

        
        function bool = iszero(X)
            bool = isempty(X.var);
        end
    end
end