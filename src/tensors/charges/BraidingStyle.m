classdef BraidingStyle
    % BraidingStyle - The braiding behaviour of charges.
    %   This represents the possibilities for the braiding of two charges.
    %
    %   Abelian - Trivial braiding with trivial twist.
    %   Bosonic - Symmetric braiding with trivial twist.
    %   Fermionic - Symmetric braiding with non-trivial twist.
    %   Anyonic - Arbitrary braiding and twist.
    %   None - No braiding defined.
    
    enumeration
        Abelian
        Bosonic
        Fermionic
        Anyonic
        None
    end
    
    methods
        function bool = issymmetric(style)
            bool = style == BraidingStyle.Abelian || style == BraidingStyle.Bosonic;
        end
        
        function c = and(a, b)
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