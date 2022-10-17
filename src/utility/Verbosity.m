classdef Verbosity < uint8
    % Verbosity level enumeration class.
    
    enumeration
        off           (0)               % No information
        warn          (1)               % Information on failure
        conv          (2)               % Convergence information
        iter          (3)               % Information about each iteration
        detail        (4)               % Detailed information about each iteration
        diagnostics   (intmax('uint8')) % all possible information
    end
end

