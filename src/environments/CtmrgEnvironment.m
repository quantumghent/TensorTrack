classdef CtmrgEnvironment
    %CTMRGENVIRONMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        corners
        edges
    end
    
    methods
        
        function obj = CtmrgEnvironment(varargin)

            if nargin == 0, return; end

            if iscell(varargin{1})
                obj.corners = varargin{1};
                obj.edges = varargin{2};
            end
            
            obj.corners = cellfun(@(x)x ./ norm(x), obj.corners,'UniformOutput',false);
            obj.edges = cellfun(@(x)x ./ norm(x), obj.edges,'UniformOutput',false);
        end
        
        function out = rot90(in)
            %rotate unit cell
            out = CtmrgEnvironment(flip(permute(in.corners,[1,3,2]),3),flip(permute(in.edges,[1,3,2]),3));
            %relabel corners and edges accordingly
            in = out;
            for i = 1:4
                out.corners(i,:,:) = in.corners(prev(i,4),:,:);
                out.edges(i,:,:) = in.edges(prev(i,4),:,:);
            end

        end

        function chi = bond_dimensions(obj)
            chi = reshape(cellfun(@(x) dims(x,1),obj.corners(1,:,:)), size(obj.corners,2:3));
        end

        function h = height(obj)
            % vertical period over which the peps is translation invariant.
            h = height(obj.corners{1,:,:});
        end
        
        function w = width(obj)
            % horizontal period over which the peps is translation invariant.
            w = width(obj.corners{1,:,:});
        end

    end
end

