classdef (InferiorClasses = {?Tensor, ?MpsTensor, ?SparseTensor}) PepsSandwich
    % Data structure representing a pair of PEPS tensors in an overlap, which behave as an
    % MPO tensor.
    
    properties
        top PepsTensor
        bot PepsTensor
    end
    
    
    methods
        function T = PepsSandwich(top, bot)
            arguments
                top
                bot = conj(top)
            end
            
            if ~isa(top, 'PepsTensor')
                top = PepsTensor(top);
            end
            if ~isa(bot, 'PepsTensor')
                bot = PepsTensor(bot);
            end
            
            T.top = top;
            T.bot = bot;
        end
        
        function T = rot90(T)
            T.top = tpermute(T.top, [1, 3, 4, 5, 2], rank(T.top));
            T.bot = tpermute(T.bot, [1, 3, 4, 5, 2], rank(T.bot));
        end
        
        function T = rot270(T)
            T.top = tpermute(T.top, [1, 5, 2, 3, 4], rank(T.top));
            T.bot = tpermute(T.bot, [1, 5, 2, 3, 4], rank(T.bot));
        end
        
        function T = transpose(T)
            % TODO: test if this is correct when inserted into FiniteMpo
            top = T.top;
            bot = T.bot;
            T.top = tpermute(bot, [1, 4, 5, 2, 3], rank(bot));
            T.bot = tpermute(top, [1, 4, 5, 2, 3], rank(top));
        end
        
        function T = ctranspose(T)
            % TODO: test if this is correct when inserted into FiniteMpo
            T.top = conj(tpermute(T.top, [1, 4, 5, 2, 3], rank(T.top)));
            T.bot = conj(tpermute(T.bot, [1, 4, 5, 2, 3], rank(T.bot)));
        end
        
        function s = pspace(T)
            s = domainspace(T)';
        end
        
        function s = domainspace(T)
            % flipped for consistency
            s = [northvspace(T.bot), northvspace(T.top)];
        end
        
        function s = codomainspace(T)
            s = [southvspace(T.top), southvspace(T.bot)];
        end
        
        function s = rightvspace(T)
            % flipped for consistency
            s = [eastvspace(T.bot), eastvspace(T.top)];
        end
        
        function s = leftvspace(T)
            s = [westvspace(T.top), westvspace(T.bot)];
        end
        
        function v = applychannel(T, L, R, v)
            arguments
                T PepsSandwich
                L MpsTensor
                R MpsTensor
                v
            end
            auxlegs_v = nspaces(v) - 4;
            auxlegs_l = nspaces(L) - 4;
            auxlegs_r = nspaces(R) - 4;
            newrank = rank(v); newrank(2) = newrank(2) + auxlegs_l + auxlegs_r;
            
            v = contract(...
                v, [1, 4, 3, 7, (-(1:auxlegs_v) - 4 - auxlegs_l)], ...
                L, [-1, 2, 5, 1, (-(1:auxlegs_l) - 4)], ...
                T.top.var, [6, 5, -2, 8, 4], ...
                T.bot.var, [6, 2, -3, 9, 3], ...
                R, [7, 8, 9, -4, (-(1:auxlegs_r) - 4 - auxlegs_l - auxlegs_v)], ...
                'Rank', newrank);
        end
    end
end

