classdef FusionStyle < uint8
    % The fusion product behaviour of charges.
    %
    %   Enumeration class that encodes the possible behavior for the decomposition of the
    %   fusion product of two charges.
    % 
    %   - :code:`Unique`: Single unique output.
    %   - :code:`Simple`: Multiple unique outputs.
    %   - :code:`Generic`: Multiple outputs.
    
    enumeration
        Unique      (0)
        Simple      (1)
        Generic     (2)
    end
    
    methods
        function bool = hasmultiplicity(style)
            % Determine whether a given fusionstyle admits fusion multiplicities.
            %
            % Returns :code:`true` for :code:`FusionStyle.Generic` and :code:`false` for all
            % other styles.
            bool = style == FusionStyle.Generic;
        end
        
        function style = and(style1, style2)
            % Determine the fusion style for a direct product of charges. This effectively
            % boils down to returning the least specific style.
            % 
            % Arguments
            % ---------
            % style1 : :class:`.FusionStyle`
            %   fusion style of first charge in direct product
            % style2 : :class:`.FusionStyle`
            %   fusion style of second charge in direct product
            %
            % Returns
            % -------
            % style : :class:`.FusionStyle`
            %   fusion style of direct product charge
            
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