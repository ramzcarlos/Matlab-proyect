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


prompt = 'Intervalo de tiempo (delta_t): \n';
% lo guarda en la variable delta_t
delta_t = input(prompt);

%carga
prompt = 'Periodo del excitacion (seg) T: \n';
% lo guarda en la variable T_e
T_e = input(prompt);
prompt = 'Valor de Fuerza: \n';
% lo guarda en la variable F
F = input(prompt);


omega=2*pi/T;
omega_e=2*pi/m;
omega_d=omega*(sqrt(1-(epsilon*epsilon)));

%Condiciones iniciales
x0=0;
x1=0;
G=0;
H=-0.253303;
A=0.253303;
B=0.012681;


T_s=0:0.01:6;
T_s=T_s(:);
tam_ts=size(T_s);
tam_ts=tam_ts(1);

M=zeros(tam_ts, 3);

for i=1:tam_ts
    for j=1:3
        if i==1
            if j==1
                M(i,j)=0;
            elseif j==2
                M(i,j)=0;
            else
                M(i,j)=((exp(-epsilon*omega*T_s(i)))*(-A*(omega_d*omega_d)*cos(omega_d*T_s(i))-B*(omega_d*omega_d)*sin(omega_d*T_s(i))))+((2*(-epsilon*omega)*(exp(-epsilon*omega*T_s(i))))*(-A*omega_d*sin(omega_d*T_s(i))+B*omega_d*cos(omega_d*T_s(i))))+(((epsilon*epsilon)*(omega*omega)*(exp(-epsilon*omega*T_s(i))))*(A*cos(omega_d*T_s(i))))*(A*cos(omega_d*T_s(i))+B*sin(omega_d*T_s(i)))-(G*(omega_e*omega_e)*sin(omega_e*T_s(i)))-(H*(omega_e*omega_e)*cos(omega_e*T_s(i)));
            end
        else
            if j==1
                M(i,j)=((exp((-epsilon)*omega*T_s(i)))*((A*cos(omega_d*T_s(i)))+(B*sin(omega_d*T_s(i))))+(G*sin(omega_e*T_s(i)))+(H*cos(omega_e*T_s(i))));
            elseif j==2
                M(i,j)=((exp((-epsilon)*omega*T_s(i)))*((-A*omega_d*sin(omega_d*T_s(i)))+(B*omega_d*cos(omega_d*T_s(i)))))-((epsilon*omega*(exp((-epsilon)*omega*T_s(i))))*((A*cos(omega_d*T_s(i)))+(B*sin(omega_d*T_s(i)))))+(G*omega_e*cos(omega_e*T_s(i)))-(H*omega_e*sin(omega_e*T_s(i)));
            else
                M(i,j)=((exp(-epsilon*omega*T_s(i)))*(-A*(omega_d*omega_d)*cos(omega_d*T_s(i))-B*(omega_d*omega_d)*sin(omega_d*T_s(i))))+((2*(-epsilon*omega)*(exp(-epsilon*omega*T_s(i))))*(-A*omega_d*sin(omega_d*T_s(i))+B*omega_d*cos(omega_d*T_s(i))))+(((epsilon*epsilon)*(omega*omega)*(exp(-epsilon*omega*T_s(i))))*(A*cos(omega_d*T_s(i))+B*sin(omega_d*T_s(i))))-(G*(omega_e*omega_e)*sin(omega_e*T_s(i)))-(H*(omega_e*omega_e)*cos(omega_e*T_s(i)));
            end
        end
    end
end                


%g%guarda archivos
%------------------------------------------------------------------

file_path_salve=pwd+"/Solucion_exacta.txt";
fileID = fopen(file_path_salve,'w');
fprintf(fileID,"id\tDesplazamiento\tVelocidad\tAceleracion\n");
formatSpec = '%f\t%f\t%f\t%f';

for i =1:tam_ts
    fprintf(fileID, formatSpec, T_s(i));
   for j =1:3
       fprintf(fileID,formatSpec,M(i,j));
   end
   fprintf(fileID,"\n");
end
fprintf(fileID,"\n");
fclose(fileID);




%------------------------------------------------------------------
%graficas

figure
plot(T_s,M(1:tam_ts, 1))
title("historia de desplazamiento")

figure
plot(T_s,M(1:tam_ts, 2))
title("historia de velocidades")

figure
plot(T_s,M(1:tam_ts, 3))
title("historia de aceleracion")





