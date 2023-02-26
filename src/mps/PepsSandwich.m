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
            top = T.top;
            bot = T.bot;
            T.top = tpermute(bot, [1, 4, 5, 2, 3], rank(bot));
            T.bot = tpermute(top, [1, 4, 5, 2, 3], rank(top));
        end
        
        function T = ctranspose(T)
            T.top = tpermute(conj(T.top), [1, 2, 5, 4, 3], rank(T.top));
            T.bot = tpermute(conj(T.bot), [1, 2, 5, 4, 3], rank(T.bot));
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
        
        function v = applympo(varargin)
            assert(nargin >= 3)
            v = varargin{end};
            R = varargin{end-1};
            L = varargin{end-2};
            O = varargin(1:end-3);
            W = length(O);
            
            auxlegs_v = nspaces(v) - (2*W+2);
            if isa(L, 'MpsTensor')
                auxlegs_l = L.alegs;
            else
                auxlegs_l = 0;
            end
            if isa(R, 'MpsTensor')
                auxlegs_r = R.alegs;
            else
                auxlegs_r = 0;
            end
            auxlegs_extra = auxlegs_l + auxlegs_r;
            
            inds_top = arrayfun(@(x) [  5 + 5*(x-1), ...
                                        2 + 5*(x-1), ...
                                        -(2 + 2*(x-1)), ...
                                        2 + 5*x, ...
                                        3 + 5*(x-1)], ...
                                        1:W,  'UniformOutput', false);
            inds_bot = arrayfun(@(x) [  5 + 5*(x-1), ...
                                        4 + 5*(x-1), ...
                                        -(3 + 2*(x-1)), ...
                                        4 + 5*x, ...
                                        6 + 5*(x-1)], ...
                                        1:W,  'UniformOutput', false);
            tops = cellfun(@(x) x.top.var, O, 'UniformOutput', false);
            bots = cellfun(@(x) x.bot.var, O, 'UniformOutput', false);
            Oargs = [tops; inds_top; bots; inds_bot];
            v = contract(v, [1, reshape([3, 6]' + 5*(0:W-1), 1, []), 3 + 5*W, (-(1:auxlegs_v) - (2*W+2) - auxlegs_l)], ...
                L, [-1, 4, 2, 1, (-(1:auxlegs_l) - (2*W+2))], ...
                Oargs{:}, ...
                R, [3 + 5*W, 2 + 5*W, 4 + 5*W, -(2*W + 2), (-(1:auxlegs_r) - (2*W+2) - auxlegs_l - auxlegs_v)], ...
                'Rank', rank(v) + [0 auxlegs_extra]);
        end
        
        function t = MpoTensor(T)
            fuse_east = Tensor.eye(prod(leftvspace(T)), leftvspace(T));
            fuse_south = Tensor.eye(prod(codomainspace(T)), codomainspace(T));
            fuse_west = Tensor.eye(prod(rightvspace(T))', rightvspace(T));
            fuse_north = Tensor.eye(prod(domainspace(T))', domainspace(T));
            t = MpoTensor(contract(...
                T.top.var, [1, 2, 4, 6, 8], ...
                T.bot.var, [1, 3, 5, 7, 9], ...
                fuse_east, [-1, 3, 2], ...
                fuse_south, [-2, 5, 4], ...
                fuse_west, [-3, 6, 7], ...
                fuse_north, [-4, 8, 9], ...
                'Rank', [2, 2]));
             % TODO: test properly
        end
    end
end

