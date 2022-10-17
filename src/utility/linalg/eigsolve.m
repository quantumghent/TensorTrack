function [V, D, flag] = eigsolve(A, v, howmany, which, kwargs)
%EIGSOLVE Summary of this function goes here
%   Detailed explanation goes here

arguments
    A
    v
    howmany = 1
    which {mustBeMember(which, {'lm', 'largestabs', 'lr', 'largestreal', ...
        'li', 'largestimag', 'sm', 'smallestabs', 'sr', 'smallestreal', ...
        'si', 'smallestimag'})} = 'lm'
    kwargs.Tol = eps(underlyingType(v))^(3/4)
    kwargs.Algorithm = 'Arnoldi'
    kwargs.MaxIter = 100
    kwargs.KrylovDim {mustBeGreaterThan(kwargs.KrylovDim, howmany)} = 20
    kwargs.DeflateDim
    kwargs.ReOrth = 2
    kwargs.NoBuild
    kwargs.Verbosity = Verbosity.warn
end

% Input validations
if ~isfield(kwargs, 'NoBuild')
    kwargs.NoBuild = ceil(kwargs.KrylovDim / 10);
end
if ~isfield(kwargs, 'DeflateDim')
    kwargs.DeflateDim = max(round(3/5 * kwargs.KrylovDim), howmany);
else
    assert(kwargs.DeflateDim < kwargs.KrylovDim, 'eigsolve:argerror', ...
        'Deflate size should be smaller than krylov dimension.')
end
if ~isa(A, 'function_handle')
    A = @(x) A * x;
end
if norm(v) < eps(underlyingType(v))^(3/4)
    error('eigsolve:inputnorm', 'starting vector should not have zero norm.');
end

% Arnoldi storage for Krylov subspace basis
V = zeros(0, kwargs.KrylovDim, 'like', v);

% Arnoldi storage for Hessenberg matrix
H = zeros(kwargs.KrylovDim, kwargs.KrylovDim, underlyingType(v));  

ctr_outer = 0;
ctr_inner = 0;
flag = 0;

while ctr_outer < kwargs.MaxIter
    ctr_outer = ctr_outer + 1;
    
    while ctr_inner < kwargs.KrylovDim  % build Krylov subspace
        ctr_inner = ctr_inner + 1;
        
        V(1:length(v), ctr_inner) = v;
        v = A(v);
        H(1:ctr_inner, ctr_inner) = V(:, 1:ctr_inner)' * v;
        v = v - V * H(:, ctr_inner);
        
        % reorthogonalize new vector
        if ctr_inner >= kwargs.ReOrth
            c = V(:, 1:ctr_inner)' * v;
            H(1:ctr_inner, ctr_inner) = H(1:ctr_inner, ctr_inner) + c;
            v = v - V(:, 1:ctr_inner) * c;
        end
        
        % normalize
        beta = norm(v, 'fro');
        v = v / beta;
        
        if ctr_inner == kwargs.KrylovDim, break; end
        
        if ctr_inner >= howmany
            invariantsubspace = beta < eps(underlyingType(beta))^(3/4);
            if invariantsubspace || ctr_inner == kwargs.KrylovDim
                break;
            end
            
            % check for convergence during subspace build
            if ~mod(ctr_inner, kwargs.NoBuild)
                [U, lambda] = eig(H(1:ctr_inner, 1:ctr_inner), 'vector');
                select = selecteigvals(lambda, howmany, which);
                conv = max(abs(beta * U(ctr_inner, select)));
                
                if conv < kwargs.Tol
                    V = V(:, 1:ctr_inner) * U(:, select);
                    D = diag(lambda(select));
                    if kwargs.Verbosity >= Verbosity.conv
                        fprintf('Conv %2d (%2d/%2d): error = %.5e.\n', ctr_outer, ...
                            ctr_inner, kwargs.KrylovDim, conv);
                    end
                    return
                end
                
                if kwargs.Verbosity >= Verbosity.detail
                    fprintf('Iter %2d (%2d/%2d):\tlambda = %.5e + %.5ei;\terror = %.5e\n', ...
                        ctr_outer, ctr_inner, kwargs.KrylovDim, ...
                        real(lambda(1)), imag(lambda(1)), conv);
                end
            end
        end
        
        H(ctr_inner + 1, ctr_inner) = beta;
    end
    
    % stopping criterium reached - irrespective of convergence
    if ctr_outer == kwargs.MaxIter || ctr_inner ~= kwargs.KrylovDim
        [U, lambda] = eig(H(1:ctr_inner, 1:ctr_inner), 'vector');
        select = selecteigvals(lambda, howmany, which);
        conv = max(abs(beta * U(ctr_inner, select)));
        V = V(:, 1:ctr_inner) * U(:, select);
        D = diag(lambda(select));
        
        if conv > kwargs.Tol
            if invariantsubspace
                fprintf('Found invariant subspace.\n');
                flag = 1;
            else
                fprintf('Reached maxiter without convergence.\n');
                flag = 2;
            end
        end
        return
    end
    
    % deflate Krylov subspace
    [U1, T] = schur(H, 'real');
    E = ordeig(T);
    select1 = false(size(E));
    select1(selecteigvals(E, kwargs.DeflateDim, which)) = true;
    [U1, T] = ordschur(U1, T, select1);
    
    V = V * U1;
    [U, lambda] = eig(T(1:kwargs.DeflateDim, 1:kwargs.DeflateDim), 'vector');
    select = selecteigvals(lambda, howmany, which);
    conv = max(abs(beta * U1(kwargs.KrylovDim, 1:kwargs.DeflateDim) * U(:, select)));
    
    % check for convergence
    if conv < kwargs.Tol
        V = V(:, 1:kwargs.DeflateDim) * U(:, select);
        D = diag(lambda(select));
        if kwargs.Verbosity >= Verbosity.conv
            fprintf('Conv %2d: error = %.5e.\n', ctr_outer, conv);
        end
        return
    end
    
    if kwargs.Verbosity >= Verbosity.iter
        fprintf('Iter %2d:\tlambda = %.5e + %.5ei;\terror = %.5e\n', ...
            ctr_outer, real(lambda(1)), imag(lambda(1)), conv);
    end
    
    % deflate Krylov subspace
    H = zeros(kwargs.KrylovDim, kwargs.KrylovDim, underlyingType(v));
    H(1:kwargs.DeflateDim, 1:kwargs.DeflateDim) = ...
        T(1:kwargs.DeflateDim, 1:kwargs.DeflateDim);
    H(kwargs.DeflateDim + 1, 1:kwargs.DeflateDim) = ...
        beta * U1(kwargs.KrylovDim, 1:kwargs.DeflateDim);
%     V(:, kwargs.DeflateDim + 1:end) = 0 * V(:, kwargs.DeflateDim + 1:end);
    ctr_inner = kwargs.DeflateDim;
end

end

function select = selecteigvals(lambda, howmany, which)

switch which
    case {'largestabs', 'lm'}
        [~, p] = sort(abs(lambda), 'descend');
    case {'largestreal', 'lr'}
        [~, p] = sort(real(lambda), 'descend');
    case {'largestimag', 'li'}
        [~, p] = sort(imag(lambda), 'descend');
    case {'smallestabs', 'sm'}
        [~, p] = sort(abs(lambda), 'ascend');
    case {'smallestreal', 'sr'}
        [~, p] = sort(real(lambda), 'ascend');
    case {'smallestimag', 'si'}
        [~, p] = sort(imag(lambda), 'ascend');
end

select = p(1:howmany);

end
