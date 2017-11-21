function [rowMap,colMap] = get_node_mapping_matrix(tri)
I = ones(1,12);
colMap = [I*tri(1),I*tri(2),I*tri(3),I*tri(4),I*tri(5),I*tri(6),I*tri(7),I*tri(8),I*tri(9),I*tri(10),I*tri(11),I*tri(12),];
rowMap = [tri',tri',tri',tri',tri',tri',tri',tri',tri',tri',tri',tri'];
end
