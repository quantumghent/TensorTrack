classdef FusionStyle < uint8
    % FusionStyle - The fusion product behaviour of charges.
    %   This represents the possibilities for the decomposition of the fusion product of two charges.
    %
    %   Unique - Single unique output.
    %   Simple - Multiple unique outputs.
    %   Generic - Multiple outputs.
    
    enumeration
        Unique      (0)
        Simple      (1)
        Generic     (2)
    end
    
    methods
        function bool = hasmultiplicity(style)
            bool = style == FusionStyle.Generic;
        end
        
        function style = and(style1, style2)
            % Determine the fusionstyle for a direct product of charges. This effectively
            % boils down to returning the least specific style.
            
            if style1 == FusionStyle.Generic || style2 == FusionStyle.Generic
                style = FusionStyle.Generic;
                return
            end
            
            if style1 == FusionStyle.Simple || style2 == FusionStyle.Simple
                style = FusionStyle.Simple;
                return
            end
            
            style = FusionStyle.Unique;
        end
    end
end