function [A,F,inner_vertices] = remove_boundary(A,F,edge)
sort_edge = sort(unique(edge(:,1)));
j = 1; i = 1; cnt = 1;
inner_vertices = zeros(length(A)-length(edge),1);
while i <= length(A)
    if i == sort_edge(j)
        i = i+1;
        j = j+1;
    else 
        inner_vertices(cnt) = i;
        i = i+1;
        cnt = cnt+1;      
    end
end

A = A(inner_vertices,inner_vertices);
F = F(inner_vertices);