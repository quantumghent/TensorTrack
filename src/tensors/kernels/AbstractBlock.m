classdef (Abstract) AbstractBlock
    % Abstract structure for storing tensor data.
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
            % fun : function_handle
            %   initialising function for the tensor data, with signature
            %   `data = fun(dims)` where dims is a row vector of dimensions.
            %
            % codomain : AbstractSpace
            %   vector of vector spaces that form the codomain.
            %
            % domain : AbstractSpace
            %   vector of vector spaces that form the domain.
            %
            % Returns
            % -------
            % X : AbstractBlock
            %   tensor data.
            
            if isa(codomain, 'CartesianSpace') || isa(codomain, 'ComplexSpace') || ...
                    isa(domain, 'CartesianSpace') || isa(domain, 'ComplexSpace')
                X = TrivialBlock(codomain, domain);
                return
            end
            
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
%             if ~strcmp(class(charges(codomain)), class(med(1).charge))
%                 bla
%             end
            X = med;
        end
    end
    
    
    %% Required methods
    methods
        function Y = axpby(a, X, b, Y, p, map)
            % (Abstract) Compute ```Y = permute(X, p) .* a + Y .* b```.
            % This method is the computationally critical method of this class, thus has
            % special cases for scalar multiplication (a == 0), addition (nargin == 4), and
            % various optimizations when a == 1, b == 0 | b == 1. Additionally, this method
            % should not be called directly as it should not perform any error checks.
            %
            % Arguments
            % ---------
            % a : double
            %   scalar to multiply with X.
            %
            % X : :class:`AbstractBlock`
            %   list of source blocks.
            %
            % b : double
            %   scalar to multiply with Y.
            %
            % Y : :class:`AbstractBlock`
            %   list of destination blocks.
            %
            % p : int
            %   permutation vector for X.
            %
            % map : (sparse) double
            %   coefficient matrix for permuting X.
            %
            % Returns
            % -------
            % Y : :class:`AbstractBlock`
            %   Result of computing Y = permute(X, p) .* a + Y .* b.
            
            error('This method should be overloaded.');
        end
        
        function [mblocks, mcharges] = matrixblocks(b)
            % (Abstract) Extract a list of coupled matrix blocks.
            %
            % Arguments
            % ---------
            % b : AbstractBlock
            %   list of input data.
            %
            % Returns
            % -------
            % mblocks : cell
            %   list of non-zero coupled matrix blocks, sorted according to its charge.
            %
            % mcharges : AbstractCharge
            %   list of coupled charges.
            
            error('This method should be overloaded.');
        end
        
        function C = mul(C, A, B, a, b)
            % (Abstract) Compute ```C = (A .* a) * (B .* b)```.
            % Compute the matrix product of two source tensor structures, and store the
            % result in the destination tensor structure. This method should not perform any
            % error checks.
            %
            % Arguments
            % ---------
            % C : AbstractBlock
            %   location to store the result
            %
            % A : AbstractBlock
            %   first matrix factor
            %
            % B : AbstractBlock
            %   second matrix factor
            %
            % a : double = 1
            %   first scalar factor
            %
            % b : double = 1
            %   second scalar factor
            %
            % Returns
            % -------
            % C : AbstractBlock
            %   Result of computing C = (A .* a) * (B .* b)
            
            error('This method should be overloaded.');
        end
        
        function tblocks = tensorblocks(b)
            % (Abstract) Extract a list of uncoupled tensor blocks.
            %
            % Arguments
            % ---------
            % b : AbstractBlock
            %   list of input data.
            %
            % Returns
            % -------
            % tblocks : cell
            %   list of non-zero uncoupled tensor blocks, sorted according to the coupled
            %   charge and then in column-major order according to the uncoupled charges.
            
            error('This method should be overloaded.');
        end
    end
    
    
    %% Optional methods
    methods
        function Y = axpy(a, X, Y, p, map)
            % Compute ```Y = permute(X, p) .* a + Y```.
            % This method is a convenience method that automatically falls back on
            % ```Y = axpby(a, X, 1, Y, p, map)```, but can be overloaded if the additional
            % efficiency is desired.
            %
            % Arguments
            % ---------
            % a : double
            %   scalar to multiply with X.
            %
            % X : AbstractBlock
            %   list of source blocks.
            %
            % Y : AbstractBlock
            %   list of destination blocks.
            %
            % p : int
            %   permutation vector for X.
            %
            % map : (sparse) double
            %   coefficient matrix for permuting X.
            %
            % Returns
            % -------
            % Y : AbstractBlock
            %   Result of computing Y = permute(X, p) .* a + Y.
            
            Y = axpby(a, X, 1, Y, p, map);
        end
        
        function Y = minus(X, Y)
            % Subtraction of X and Y.
            %
            % Usage
            % -----
            % Y = minus(X, Y)
            % Y = X - Y
            %
            % Arguments
            % ---------
            % X : AbstractBlock
            %   first list of input data.
            %
            % Y : AbstractBlock
            %   second list of input data.
            %
            % Returns
            % -------
            % Y : AbstractBlock
            %   list of output matrices.
            
            Y = axpby(1, X, -1, Y);
        end
        
        function Y = plus(X, Y)
            % Addition of X and Y.
            %
            % Usage
            % -----
            % Y = plus(X, Y)
            % Y = X + Y
            %
            % Arguments
            % ---------
            % X : AbstractBlock
            %   first list of input data.
            %
            % Y : AbstractBlock
            %   second list of input data.
            %
            % Returns
            % -------
            % Y : MatrixBlock
            %   list of output data.
            
            Y = axpby(1, X, 1, Y);
        end
        
        function Y = times(Y, a)
            % Scalar multiplication of Y and a.
            %
            % Usage
            % -----
            % Y = times(Y, a)
            % Y = Y .* a
            %
            % Arguments
            % ---------
            % Y : AbstractBlock
            %   list of input data.
            %
            % a : double
            %   scalar factor.
            %
            % Returns
            % -------
            % Y : AbstractBlock
            %   list of output data.
            
            Y = axpby(0, [], a, Y);
        end
        
        function Y = rdivide(Y, a)
            % Scalar division of Y and a.
            %
            % Usage
            % -----
            % Y = rdivide(Y, a)
            % Y = Y ./ a
            %
            % Arguments
            % ---------
            % Y : AbstractBlock
            %   list of input data.
            %
            % a : double
            %   scalar factor.
            %
            % Returns
            % -------
            % Y : AbstractBlock
            %   list of output data.
            
            Y = axpby(0, [], 1/a, Y);
        end
        
        function A = uplus(A)
            % Unary plus. Equivalent to making a copy.
            %
            % Usage
            % -----
            % A = uplus(A)
            % A = +A
            %
            % Arguments
            % ---------
            % A : AbstractBlock
            %   list of input data.
            %
            % Returns
            % -------
            % A : AbstractBlock
            %   list of output data.
            
        end
        
        function Y = uminus(Y)
            % Unary minus. Computes the additive inverse.
            %
            % Usage
            % -----
            % A = uminus(A)
            % A = -A
            %
            % Arguments
            % ---------
            % A : MatrixBlock
            %   list of input matrices.
            %
            % Returns
            % -------
            % A : MatrixBlock
            %   list of output matrices.
            
            Y = Y .* (-1);
        end
        
        function t = underlyingType(X)
            % The scalar type of the tensor.
            %
            % Arguments
            % ---------
            % X : :class:`AbstractBlock`
            %   input block.
            % 
            % Returns
            % -------
            % type : char
            %   the scalar type of the data, which is one of the following:
            %
            %   - 'single' or 'double'
            %   - 'logical'
            %   - 'int8', 'int16', 'int32' or 'int64'
            %   - 'uint8', 'uint16', 'uint32' or 'uint64'
            
            error('tensors:AbstractMethod', 'This method should be overloaded.');
        end
    end
end