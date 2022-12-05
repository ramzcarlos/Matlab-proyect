clear;
clc;

%***************************************************
%solicita al usuario los valores necesarios para realizar los calculos
prompt = 'Valor de la masa (m): \n';
% lo guarda en la variable m
m = input(prompt);

prompt = 'Periodo del oscilador (seg) T: \n';
% lo guarda en la variable T
T = input(prompt);

prompt = 'Valor del amortiguamiento (epsilon): \n';
% lo guarda en la variable epsilon
epsilon = input(prompt);
epsilon=epsilon/100;

prompt = 'Valor de la amplitud (F): \n';
% lo guarda en la variable F
F = input(prompt);

prompt = 'Intervalo de tiempo (delta_t): \n';
% lo guarda en la variable delta_t
delta_t = input(prompt);

%-------------------------------------------------------------------------
% comienzan a realizar los caluculos

beta=1/4;
gama=1/2;
%Rigidez
k=(((2*pi)/T)*((2*pi)/T))*m;
c=2*epsilon*(sqrt(k*m));
omega=sqrt(k/m);

a1=((1/(beta*delta_t*delta_t))*m)+((gama/(beta*delta_t))*c);
a2=(((1/(beta*delta_t))*m)+(((gama/beta)-1)*c));
a3=(((1/(2*beta))-1)*m)+((delta_t*((gama/(2*beta))-1)*c));

k1=k+a1;
p0=0;
u=0;
u1=0;

T_s=0:0.01:6;
T_s=T_s(:);
tam_ts=size(T_s);
tam_ts=tam_ts(1);

M=zeros(tam_ts, 5);

for i=1:tam_ts
    for j=1:5
        if i==1
            if j==1
                M(i,j)=(F/m)*sin(omega*T_s(i));
            elseif j==2
                M(i,j)=T_s(i);
            elseif j==3
                M(i,j)=u;
            elseif j==4
                M(i,j)=u1;
            else
                M(i,j)=(M(i,1)-(c*M(i,3)-(k*M(i,2))))/m;
            end
        else
            if j==1
                M(i,j)=(F/m)*sin(omega*T_s(i));
            elseif j==2
                M(i,j)=M(i,1)+(a1*M(i-1,3))+(a2*M(i-1,4))+(a3*M(i-1,5));
            elseif j==3
                M(i,j)=M(i,2)/k1;
            elseif j==4
                M(i,j)=(2/delta_t)*(M(i,3)-M(i-1,3))-M(i-1,4);
            else
                M(i,j)=((4/(delta_t*delta_t))*(M(i,3)-M(i-1,3)))-((4/delta_t)*(M(i-1,4)))-M(i-1,5);
            end
        end

    end
end

%g%guarda archivos
%------------------------------------------------------------------

file_path_salve=pwd+"/Aceleracion_contante_rigidez.txt";
fileID = fopen(file_path_salve,'w');
fprintf(fileID,"id\tDesplazamiento\tVelocidad\tAceleracion\n");
formatSpec = '%f\t%f\t%f\t%f';

for i =1:tam_ts
    fprintf(fileID, formatSpec, T_s(i));
   for j =3:5
       fprintf(fileID,formatSpec,M(i,j));
   end
   fprintf(fileID,"\n");
end
fprintf(fileID,"\n");
fclose(fileID);




%------------------------------------------------------------------
%graficas

figure
plot(T_s,M(1:tam_ts, 3))
title("historia de desplazamiento")

figure
plot(T_s,M(1:tam_ts, 4))
title("historia de velocidades")

figure
plot(T_s,M(1:tam_ts, 5))
title("historia de aceleracion")



