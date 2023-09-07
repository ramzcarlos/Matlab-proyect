clear;
clc;
prompt= 'Ingrese el Número de filas de refuerzo a calcular: \n';
n=input(prompt);

v=zeros([n 1]);
v=v(:);

dv=zeros([n 1]);
dv=v(:);
dd=zeros([n 1]);
dd=v(:);


for i=1:n
    msg='Ingrese el Número de varillas del nivel  ';
    msg1=int2str(i);
    msg2='-';
    msg3=': \n';
    msg=strcat(msg, msg2, msg1, msg3);
    prompt= msg;
    dv(i)=input(prompt);
    msg='Ingrese el diametro de varillas del nivel  ';
    msg1=int2str(i);
    msg2=' ';
    msg3=' en pulgadas: \n';
    msg=strcat(msg, msg2, msg1, msg3);
    prompt= msg;
    dd(i)=input(prompt);
end