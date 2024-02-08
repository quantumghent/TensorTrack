classdef BraidingStyle
    % The braiding behaviour of charges.
    %
    %   Enumeration class that encodes the possibilities for the braiding of two charges:
    % 
    %   - :code:`Abelian`: Trivial braiding with trivial twist.
    %   - :code:`Bosonic`: Symmetric braiding with trivial twist.
    %   - :code:`Fermionic`: Symmetric braiding with non-trivial twist.
    %   - :code:`Anyonic`: Arbitrary braiding and twist.
    %   - :code:`None`: No braiding defined.
    
    enumeration
        Abelian
        Bosonic
        Fermionic
        Anyonic
        None
    end
    
    methods
        function bool = istwistless(style)
            % Determine whether a given braiding style has a trivial twist.
            %
            % Returns :code:`true` for :code:`BraidingStyle.Abelian` and
            % :code:`BraidingStyle.Bosonic` and :code:`false` for all other styles.
            bool = style == BraidingStyle.Abelian || style == BraidingStyle.Bosonic;
        end
        
        function bool = issymmetric(style)
            % Determine whether a given braiding style is symmetric.
            %
            % Returns :code:`true` for all styles except :code:`BraidingStyle.Anyonic` and
            % :code:`BraidingStyle.None`.
            bool = style == BraidingStyle.Abelian || style == BraidingStyle.Bosonic || ...
                style == BraidingStyle.Fermionic;
        end
        
        function c = and(a, b)
            % Determine the braiding style for a direct product of charges. This effectively
            % boils down to returning the least specific style.
            % 
            % Arguments
            % ---------
            % style1 : :class:`.BraidingStyle`
            %   fusion style of first charge in direct product
            % style2 : :class:`.BraidingStyle`
            %   fusion style of second charge in direct product
            %
            % Returns
            % -------
            % style : :class:`.BraidingStyle`
            %   fusion style of direct product charge
            if any([a b] == BraidingStyle.None)
                c = BraidingStyle.None;
                
            elseif any([a b] == BraidingStyle.Anyonic)
                c = BraidingStyle.Anyonic;
                
            elseif any([a b] == BraidingStyle.Fermionic)
                c = BraidingStyle.Fermionic;
                
            elseif any([a b] == BraidingStyle.Bosonic)
                c = BraidingStyle.Bosonic;
                
            else
                c = BraidingStyle.Abelian;
                
            end
        end
    end
end