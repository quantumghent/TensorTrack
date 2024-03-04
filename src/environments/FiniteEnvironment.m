classdef FiniteEnvironment
    % Data structure for managing the environments in finite MPS algorithms.
    %
    % Properties
    % ----------
    % GL : :class:`cell` of :class:`.MpsTensor`
    %   environment tensors corresponding to the left-gauged part of the MPS.
    %
    % GR : :class:`cell` of :class:`.MpsTensor`
    %   environment tensors corresponding to the right-gauged part of the MPS.
    %
    % Todo
    % ----
    % Document
    
    properties
        GL
        GR
    end
    
    methods
        function envs = FiniteEnvironment(varargin)
            if nargin == 0, return; end
            if nargin == 2
                envs.GL = varargin{1};
                envs.GR = varargin{2};
                assert(isequal(size(envs.GL), size(envs.GR)));
                return
            end
            error('undefined syntax');
        end
        
        function envs = movegaugecenter(envs, mpo, mps1, mps2, pos)
            for i = 2:pos
                if ~isempty(envs.GL{i}), continue; end
                T = transfermatrix(mpo, mps1, mps2, i-1);
                envs.GL{i} = apply(T, envs.GL{i-1});    
            end
            for i = length(mps1):-1:(pos+1)
                if ~isempty(envs.GR{i}), continue; end
                T = transfermatrix(mpo, mps1, mps2, i).';
                envs.GR{i} = apply(T, envs.GR{i+1});
            end
        end
        
        function envs = invalidate(envs, pos)
            envs.GL(pos+1:end) = cell(1, length(envs.GL) - pos);
            envs.GR(1:pos) = cell(1, pos);
        end
    end
    
%     methods (Static)
%         function envs = initialize(mpo, mps1, mps2)
%             arguments
%                 mpo
%                 mps1
%                 mps2 = mps1
%             end
%             
%             GL = cell(1, length(mpo) + 1);
%             GR = cell(1, length(mpo) + 1);
%             
%             GL{1} = mpo.L;
%             GR{end} = mpo.R;
%             
%             envs = FiniteEnvironment(GL, GR);
%         end
%     end
end

