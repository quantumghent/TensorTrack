classdef Arrow < logical
% Enumeration class reprenting the possible directions of a tensor leg:
%
% - :code:`in`: Incoming tensor leg, encoded as a logical :code:`true`
% - :code:`out`: Outgoing tensor leg, encoded as a logical :code:`false`

enumeration
    in (true)
    out (false)
end

end
