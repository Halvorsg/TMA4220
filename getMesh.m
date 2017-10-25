function getMesh()
cnt = 1;
fileID = fopen('sphere.msh','r');
points = importdata('sphere.msh',' ',5);
str = '0';
p = zeros(1000,4);
i = 1;
while ~isequal(str,'$EndNodes')
    str = fgetl(fileID);
    A = split(str);
    if length(A) == 4
        A = str2double(A)';
        p(i,:) = A;
        i = i+1;
    end
    cnt = cnt+1;
end
Elements = zeros(ceil(length(p)*10),4);
i = 1;
cnt = 1;
while ~isequal(str,'$EndElements')
    str = fgetl(fileID);
    A = split(str);
    if length(A) == 9
        A = str2double(A)';
        Elements(i,:) = A(1:4);
        i = i+1;
    end
    cnt = cnt+1;
end
disp(i)
end