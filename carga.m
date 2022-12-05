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

prompt = 'Periodo del exitacion (seg) T: \n';
% lo guarda en la variable T
T_e = input(prompt);

prompt = 'Valor del amortiguamiento (epsilon): \n';
% lo guarda en la variable epsilon
epsilon = input(prompt);
epsilon=epsilon/100;


prompt = 'Intervalo de tiempo (delta_t): \n';
% lo guarda en la variable delta_t
delta_t = input(prompt);

prompt = 'Frecuencia (seg) F: \n';
% lo guarda en la variable F
F = input(prompt);

%Carga

omega=2*pi/T;
omega_e=2*pi/T_e;
omega_d=omega*sqrt(1-(epsilon*epsilon));

%Constates
S=0.062712;
C=0.998032;
E=0.996863;
%Constantes de posicion
A0=0.998031;
A1=0.009962;
A2=0.000050;
A3=0.000000166;
%Constantes de velocidad
B0=-0.393288;
B1=0.991771;
B2=0.009962;
B3=0.000050;

T_s=0:0.01:6;
T_s=T_s(:);
tam_ts=size(T_s);
tam_ts=tam_ts(1);

M=zeros(tam_ts, 5);

for i=1:tam_ts
    for j=1:5
        if i==1
            if j==1
                M(i,j)=F*sin(omega_e*T_s(i));
            elseif j==2
                M(i,j)=0;
            elseif j==3
                M(i,j)=0;
            elseif j==4
                M(i,j)=0;
            else
                M(i,j)=0;
            end
        else
            if j==1
                M(i,j)=F*sin(omega_e*T_s(i));
            elseif j==2
                M(i,j)=(M(i,1)-M(i-1,1))/delta_t;
            elseif j==3
                M(i,j)=M(i-1,3)*A0+A1*M(i-1,4)+A2*M(i-1,1)+A3*M(i,2);
            elseif j==4
                M(i,j)=B0*M(i-1,3)+B1*M(i-1,4)+B2*M(i-1,1)+B3*M(i,2);
            else
                M(i,j)=M(i,1)/m-2*epsilon*omega*M(i,4)-(omega_d*omega_d)*M(i,3);
            end
        end
    end
end

%g%guarda archivos
%------------------------------------------------------------------

file_path_salve=pwd+"/Carga_lineal.txt";
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




