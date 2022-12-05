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


%condiciones iniciales
omega=(2*pi)/T;

x=0;
x1=0;
x2=0;
F=1;

T_s=0:0.01:6;
T_s=T_s(:);
tam_ts=size(T_s);
tam_ts=tam_ts(1);

M=zeros(tam_ts, 4);
for i=1:tam_ts
    for j=1:4
        if i==1
            if j==1
                M(i,j)=(x2/m)*sin(omega*T_s(i));
            elseif j==2
                M(i,j)=0;
            elseif j==3
                M(i,j)=0;
            else
                M(i,j)=0;
            end
        else
            if j==1
                M(i,j)=(F/m)*sin(omega*T_s(i));
            elseif j==2
                M(i,j)=(M(i,1)-M(i-1,1)-(omega*omega)*M(i-1,3)*delta_t-M(i-1,2)*(epsilon*omega*delta_t-1+((omega*omega)*(delta_t*delta_t))/4))/(1+epsilon*omega*delta_t+((omega*omega)*(delta_t*delta_t))/4);
            elseif j==3
                M(i,j)=M(i-1,3)+(delta_t/2)*(M(i,2)+M(i-1,2));
            else
                M(i,j)=M(i-1,4)+M(i-1,3)*delta_t+((delta_t*delta_t)/4)*(M(i,2)+M(i-1,2));
            end
        end
    end
end


%g%guarda archivos
%------------------------------------------------------------------

file_path_salve=pwd+"/Aceleracion_contante_frecuencia.txt";
fileID = fopen(file_path_salve,'w');
fprintf(fileID,"id\tDesplazamiento\tVelocidad\tAceleracion\n");
formatSpec = '%f\t%f\t%f\t%f';

for i =1:tam_ts
    fprintf(fileID, formatSpec, T_s(i));
   for j =2:4
       fprintf(fileID,formatSpec,M(i,j));
   end
   fprintf(fileID,"\n");
end
fprintf(fileID,"\n");
fclose(fileID);




%------------------------------------------------------------------


%graficas

figure
plot(T_s,M(1:tam_ts, 2))
title("historia de desplazamiento")

figure
plot(T_s,M(1:tam_ts, 3))
title("historia de velocidades")

figure
plot(T_s,M(1:tam_ts, 4))
title("historia de aceleracion")