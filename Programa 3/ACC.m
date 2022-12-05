function [ACC]=ACC(sismo,z,T,m,dt)
%Maestria en Ciencias de la Ingeniría
%Sonia Valdez Mejia
%
%                    VARIABLES
%
% z= Amortiguamiento
% T= Periodo
% m= Masa
% dt= Intervalo de tiempo
%sismo= Acelerograma (Desde archivo)
%
%                     FUNCIÓN
%
%function [ACC]=ACC(sismo,z,T,m,dt)
%
% -------------------Plantilla---------------------
%load sismo.txt
% [ACC]=ACC(sismo,0.05,0.5,1,0.01);
%---------------------------------------------------


n=length(sismo);
p=-m*sismo;
%Frecuencia natural y amortiguada
w=(2*pi)/T;
wd=w*sqrt(1-z*z);
%Constantes S C y E
s=sin(wd*dt);
c=cos(wd*dt);
e=exp(-z*w*dt);
%Constantes A0, A1, A2 y A3
a0=e*((w*z*s/wd)+c);
a1=e*s/wd;
a2=-e*((s*z/(m*w*wd))+c/(m*w^2))+1/(m*w^2);
a3=(e*s*(2*z*z-1)/(wd*w*w*m))+(2*z*(e*c-1)/(m*w^3))+(dt/(m*w^2));
%Constantes B0, B1, B2 y B3
b0=-(z*z*w*w+wd*wd)*e*s/wd;
b1=(-w*z*s/wd+c)*e;
b2=e*s*((z*z/(m*wd))+(wd/(m*w^2)));
b3=(1-e*c)/(m*w^2)+e*s*(((z-2*z^3)/(m*w*wd))-(2*z*wd)/(m*w^3));
%Valores iniciales
des(1)=0; vel(1)=0; acc(1)=0; nu(1)=0; 
%Iteraciones del cálculo de la respuesta
for i=2:n;
    nu(i,1)=(p(i)-p(i-1))/dt;
    des(i,1)=a0*des(i-1)+a1*vel(i-1)+a2*p(i-1)+a3*nu(i);
    vel(i,1)=b0*des(i-1)+b1*vel(i-1)+b2*(p(i-1)+b3*nu(i));
    acc(i,1)=(p(i)/m)-(2*z*w*vel(i))-(des(i)*w^2);
end
%Respuesta del sistema
if T==0
    acc_T(n,1)=(max(abs(sismo)));
% DES=max(abs(des));
% VEL=max(abs(vel));
ACC=max(abs(acc));
% ACC_T=max(abs(acc_T));
else
 acc_T=(acc+sismo);
% DES=max(abs(des));
% VEL=max(abs(vel));
ACC=max(abs(acc));
% ACC_T=max(abs(acc_T));
end