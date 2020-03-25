function [y] = FindZero(pn)
% Finds value of in that results in zero of function

%dbstop if error

y = fsolve(@fun,[1;1]);

    function g = fun(in)
        K = [-73596.5730890981,-1477.83273742007;...
            -1477.83273742007,-101477.832737420];
        g = f(pn + K*in) - in;
    end        
end

