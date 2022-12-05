clear;
clc;

disp('Hola carlos');
filename = '';
delimiterIn = ',';
headerlinesIn = 1;
A = load(filename);



for i =1:length(A)
  for j =1:length(A(i,:))
    if (A(i,j)>= 0.5)
      disp("bueno");
    else
      disp("malo");
    end
  end
end



length(A)
length(A(1,:))
path_file=pwd
fileID = fopen('celldata1.dat','w');
formatSpec = '%d %d';
tam=size(M);
n=tam(1);
m=tam(2);
for i =1:n
  for j =1:m
    %fprintf("%d %d \n", i, j)
    fprintf(fileID, formatSpec, M(i,j));
    
    %disp(A(i,j));
  end
  fprintf(fileID,"\n");
end
fclose(fileID);

