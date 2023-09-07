%function [t,T,dUt2,MUt,MdUt,MdUt2,AMdUt2] = Newmark(a,Amort,gama,beta,delta_t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%DATOS
Amort=0.05;
beta=0.05;
gama=0.5;
m=1000;              %Masa en (kg-s^2/m)
a=load("cu19ew.txt")*9.807; %registro de aceleraciones (m/s^2)
delta_t=0.02;       %delta t de los registros de acelarción (s)
t=0:delta_t:(length(a)-1)*delta_t; %tiempo para cada registro (s)
n=t(length(a))/delta_t;
T=0.1:0.1:8;            %Periodo en (s)
Tn=length(T);
delta_t=0.02;       %delta t de los registros de acelarción (s)
for c=1:length(Amort)   %# de caso de amortiguamiento
    
for j=1:Tn  %# de punto para el espectro
%Rigidez "k"
    k=m*(2*pi/T(j))^2;
%Amortiguamiento
Ccr=2*sqrt(k*m);
C=Ccr*Amort(c);
%Rigidez modificada
K=m/(beta*delta_t^2)+gama*C/(beta*delta_t)+k;
%Fuerza equivalente "Pt"
Pt=-m*a;
Ut(1)=0;
dUt(1)=0;
dUt2(1)=0;
    
     for i=1:n %graficas de Newmark
    %Fuerza equivalente modificada "Ptm"
    Ptm(i+1)=Pt(i+1)+m*(Ut(i)/(beta*delta_t^2)+dUt(i)/(delta_t*beta)+dUt2(i)*(1/(2*beta)-1))+C*(Ut(i)*gama/(beta*delta_t)+dUt(i)*(gama/beta-1)+dUt2(i)*delta_t/2*(gama/beta-2));
    
    %Desplazamiento "Ut"
    Ut(i+1)=Ptm(i+1)/K;
    
    %Aceleración "dUt2"
    dUt2(i+1)=(Ut(i+1)-Ut(i))/(beta*delta_t^2)-dUt(i)/(beta*delta_t)-dUt2(i)*(1/(2*beta)-1);
    
    %Velocidad "dUt"
    dUt(i+1)=dUt(i)+dUt2(i)*delta_t*(1-gama)+dUt2(i+1)*gama*delta_t;
      
     end
    
 %Valores Maximos
    MdUt2(j,c)=max(abs(dUt2));
    MdUt(j,c)=max(abs(dUt));
    MUt(j,c)=max(abs(Ut));
    AMdUt2(j,c)=max(abs(dUt2+a'));
    
end
end
%end

