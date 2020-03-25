function [pn_map] = pn_map_func(pn)
%pn_map_func(pn) returns pn that matches nonlinear map indices
    %The Result is Not Rounded as this must be accomplished in code

    top = 25;
    bottom = -25;
    
    pn_map = [(pn(1) - bottom)*(100/(top-bottom));
          (pn(2) - 236)*25];

    if (pn_map(1) <= 1)
        pn_map(1) = 2;
        disp('Table Saturated')
    end
    if (pn_map(1) >= 100)
        pn_map(1) = 100;
        disp('Table Saturated')
    end
    if (pn_map(2) <= 1)
        pn_map(2) = 2;
        disp('Table Saturated')
    end
    if (pn_map(2) >= 100)
        pn_map(2) = 100;
        disp('Table Saturated')
    end
end

