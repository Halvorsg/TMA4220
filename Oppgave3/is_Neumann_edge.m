function res = is_Neumann_edge(p1,p2)
    if (p1(1) == 0 && p1(2) <= 0) && (p2(1) == 0 && p2(2) <= 0)
        res = 1;
%         fprintf('Vertical\n')
    elseif (p1(1) >= 0 && p1(2) == 0) && (p2(1) >= 0 && p2(2) == 0)
        res = 2;
%         fprintf('Horizontal\n')
    else
        res = 0;
    end
end
        