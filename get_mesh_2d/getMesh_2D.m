function [p, Elements] = getMesh_2D(file_name)
cnt = 1;
fileID = fopen(file_name,'r');
points = importdata(file_name,' ',5);
str = '0';
p = zeros(100000,3);
i = 1;
while ~isequal(str,'$EndNodes')
    str = fgetl(fileID);
    A = split(str);
    if length(A) == 4
        A = str2double(A)';
        p(i,:) = A(2:4);
        i = i+1;
    end
    cnt = cnt+1;
end
Elements = zeros(ceil(length(p)*10),3);
i = 1;
cnt = 1;
while ~isequal(str,'$EndElements')
    str = fgetl(fileID);
    A = split(str);
    if length(A) == 8
        A = str2double(A)';
        Elements(i,:) = A(6:8);
        i = i+1;
    end
    cnt = cnt+1;
end
Elements=Elements(1:i-1,:);
p=p(1:max(max(Elements)),:);

end