% PROGRAMA PARA DIBUJAR EL DIAGRAMA DE INTERACCIÓN DE UNA COLUMNA
% COLUMNAS RECTANGULARES DE CONCRETO ARMADO
clc;
clear all;
disp('INGRESO DE DATOS')
fc=input ('Resistencia del Concreto Kg/cm2:');
fy=input ('Resistencia del Acero Kg/cm2:');
b=input ('Base de Columna cm:');
h=input ('Altura de Columna cm:');
ri=input ('Recubrimiento inferior cm:');
rs=input ('Recubrimiento superior cm:');
Eu=0.003;
Es=2000000;
Ey=fy/Es;
d=h-ri;
prompt= 'Número de niveles de refuerzo: ';
n=input(prompt);
v=zeros([n 1]);
v=v(:);
nv=zeros([n 1]);
nv=v(:);
dv=zeros([n 1]);
dv=v(:);
Cb=(Eu*d)/(Ey+Eu);%eje neutro balanceado
cen=h/2;%centroide de la sección
%betha

phi_c=0.65; % factor de compresion
phi_t=0.90; % factor de tension

if 170<fc<=280
    betha=0.85
else 
    if 280<fc<550
        betha=0.85-((0.05*(fc-280))/7);
    else
        betha=0.65
    end
end
    
for i=1:n
    msg='\nIngrese el Número de varillas del nivel  ';
    msg1=int2str(i);
    msg2='  ';
    msg3=':';
    msg=strcat(msg, msg2, msg1, msg3);
    prompt= msg;
    nv(i)=input(prompt);
    msg='Ingrese el diametro de varillas del nivel  ';
    msg1=int2str(i);
    msg2='  ';
    msg3=' en pulgadas:';
    msg=strcat(msg, msg2, msg1, msg3);
    prompt= msg;
    dv(i)=input(prompt);
    As(i)=nv(i).*(dv(i).*2.54)^2*(pi/4);
    espacio=(h-ri-rs)/(n-1);
    dis(i)=(espacio*(i-1))+rs;
    d1(i)=Cb-dis(i);
    esi(i)=(Eu.*d1(i))/Cb;
    if esi(i)<Ey
        esi(i)=esi(i);
    else
        esi(i)=Ey;
    end
    fsi(i)=esi(i)*Es;
    F(i)=fsi(i).*As(i);
    brazo(i)=(h/2)-dis(i);
    Ms(i)=F(i).*brazo(i);
    Fc=0.85*fc*betha*Cb*b;
    bCc=(h/2)-((betha*Cb)/2);
    Cc=Fc;
    M1=sum(Ms);
end
%Carga Axial Pura
Ar=sum(As);
P0=(0.85*fc*(b*h)+Ar*fy);
M0=0;
e0=M0/P0;
%Falla balanceada
As(:);
dis(:);
d1(:);
esi(:);
fsi(:);
F(:);
brazo(:);
Ms(:);
Pb=sum(F)+Cc;
Mb=M1+(Cc*bCc);
d2=zeros([n 1]);
d2=d2(:);
esi2=esi;
esi1=esi;
brazo1=brazo;
Ms2=zeros([n 1]);
Ms2=Ms2(:);
cont=1;
rM=[];
rP=[];
phirp=[];
phirm=[];
bandera=0;
for i=0.02:0.02:1
    %cosiente es c/d
    %if cont == 1
        %rM(cont)=M0;
        %rP(cont)=P0;
        %cont=cont+1;
    %end
    cosiente=i;
    c=cosiente*d;
    for j=1:n
        d2(j)=c-dis(j);
        esi1(j)=(Eu*d2(j))/c;
        if abs(esi1(j))<Ey
            esi2(j)=abs(esi1(j));
            

        else
            esi2(j)=Ey;
        end
        
        
        

        %Realiza los calculos de fsi
        fsi(j)=esi2(j)*Es;
        if esi1(j)<0
            F(j)=-(fsi(j).*As(j));
        else
            F(j)=(fsi(j).*As(j));
        end
        brazo1(j)=(h/2)-dis(j);
        Ms2(j)=F(j)*brazo1(j);
        

    end
    m2=sum(Ms2);
    Cc2=0.85*fc*betha*c*b;
    bCc2=(h/2)-((betha*c)/2);
    if i >= 0.02 && i <0.58
        disp(i)
        P=sum(F)+Cc2;
        m=m2+(Cc2*bCc2);
        rM(cont)=m;
        rP(cont)=P;
        if abs(esi1(n))<=Ey
            m=m*phi_c;
            phirm(cont)=m;
            P=P*phi_c;
            phirp(cont)=P;
        elseif abs(esi1(n))>Ey && abs(esi1(n))<(Ey+0.003)
            phi_trans= 0.656+0.25*((abs(esi1(n))-Ey)/0.003);
            m=m*phi_trans;
            phirm(cont)=m;
            P=P*phi_trans;
            phirp(cont)=P;
        else
            m=m*phi_t;
            phirm(cont)=m;
            P=P*phi_t;
            phirp(cont)=P;
        end



        
        cont=cont+1;
    elseif i==0.58
        P=sum(F)+Cc2;
        m=m2+(Cc2*bCc2);
        rM(cont)=m;
        rP(cont)=P;
        if abs(esi1(n))<=Ey
            m=m*phi_c;
            phirm(cont)=m;
            P=P*phi_c;
            phirp(cont)=P;
        elseif abs(esi1(n))>Ey && abs(esi1(n))<(Ey+0.003)
            phi_trans= 0.656+0.25*((abs(esi1(n))-Ey)/0.003);
            m=m*phi_trans;
            phirm(cont)=m;
            P=P*phi_trans;
            phirp(cont)=P;
        else
            m=m*phi_t;
            phirm(cont)=m;
            P=P*phi_t;
            phirp(cont)=P;
        end
        cont=cont+1;


        rM(cont)=Mb;
        rP(cont)=Pb;
        cont=cont+1;
    else
        P=sum(F)+Cc2;
        m=m2+(Cc2*bCc2);
        rM(cont)=m;
        rP(cont)=P;
        if abs(esi1(n))<=Ey
            m=m*phi_c;
            phirm(cont)=m;
            P=P*phi_c;
            phirp(cont)=P;
        elseif abs(esi1(n))>Ey && abs(esi1(n))<(Ey+0.003)
            phi_trans= 0.656+0.25*((abs(esi1(n))-Ey)/0.003);
            m=m*phi_trans;
            phirm(cont)=m;
            P=P*phi_trans;
            phirp(cont)=P;
        else
            m=m*phi_t;
            phirm(cont)=m;
            P=P*phi_t;
            phirp(cont)=P;
        end
        cont=cont+1;
    end


end
rM(cont)=M0*phi_t;
rP(cont)=P0*phi_t;
cont=cont+1;

rM=rM(:);
rP=rP(:);
phirm=phirm(:);
phirp=phirp(:);

plot(rM, rP,"-ok");
%plot(phirm, phirp, "--k");










%HASTA AQUÍ
Ass=nvs*(dvs*2.54)^2*pi/4;%Área de varillas del nivel 1
Asm=nvm*(dvm*2.54)^2*pi/4;%Área de varillas del nivel 2
Asi=nvi*(dvi*2.54)^2*pi/4;%Área de varillas del nivel 3
Ast=Ass+Asm+Asi;


% SE CALCULAS LOS MOMENTOS Y CARGAS NOMINALES
disp('CÁLCULOS')
disp('1) Fuerza Neta en Ton')
P0=(0.85*fc*(b*h)+Ast*fy)/1000
M0=0
e0=M0/P0
disp('2) Falla Balanceada')
d=h-rs;
Cb=(Eu*d)/(Ey+Eu)
a=0.85*Cb
% Condicionante para ver que aceros están en tensión y compresión
if Cb>h/2
    disp('2.1) Acero en tensión')
    As1=Asi; % Área en tensión
    Es1=Ey;
    fs1=Es*Es1;
    disp('2.1) Acero en Compresión')
    As2=Asm; % Área en Compresión
    Es2=Eu*(Cb-h/2)/Cb;
    fs2=Es*Es2;
    if fs2<fy
        fs2;
    else
        fs2=fy;
    end
    As3=Ass; % Área en Compresión
    Es3=Eu*(Cb-h+d)/Cb;
    fs3=Es*Es3;
    if fs3<fy
        fs3;
    else
        fs3=fy;
    end
    Pb=(0.85*fc*a*b+As2*fs2+As3*fs3-As1*fs1)/1000
    Mb=(0.85*fc*a*b*(h/2-a/2)+As3*fs3*(d-h/2)+As1*fs1*(d-h/2))/1000
    eb=Mb/Pb      
else
    disp('2.1) Acero en Tracción')
    As1=Asi; % Área en Tracción
    Es1=Eu*(d-Cb)/Cb;
    fs1=Es*Es1;
    if fs1<fy
        fs1;
    else
        fs1=fy;
    end
    As2=Asm; % Área en Tracción
    Es2=Eu*(h/2-Cb)/Cb;
    fs2=Es*Es2;
    if fs2<fy
        fs2;
    else
        fs2=fy;
    end    
    disp('2.1) Acero en Compresión')
    As3=Ass; % Área en Compresión
    Es3=Eu*(Cb-h+d)/Cb;
    fs3=Es*Es3;
    if fs3<fy
        fs3;
    else
        fs3=fy;
    end  
    Pb=(0.85*fc*a*b+As3*fs3-As1*fs1-As2*fs2)/1000
    Mb=(0.85*fc*a*b*(h/2-a/2)+As3*fs3*(d-h/2)+As1*fs1*(d-h/2))/1000
    eb=Mb/Pb
end
for i=1:1:49
    C(i)=h-d+(2*d-h)*i/50;
    a=0.85*C(i);
    if C(i)>h/2
    disp('2.1) Acero en Tracción')  
    As1=Asi; % Área en Tracción
    Es1=Eu*(d-C(i))/C(i);
    fs1=Es*Es1;
    if fs1<fy
        fs1;
    else
        fs1=fy;
    end
    disp('2.1) Acero en Compresión')
    As2=Asm; % Área en Compresión
    Es2=Eu*(C(i)-h/2)/C(i);
    fs2=Es*Es2;
    if fs2<fy
        fs2;
    else
        fs2=fy;
    end
    As3=Ass; % Área en Compresión
    Es3=Eu*(C(i)-h+d)/C(i);
    fs3=Es*Es3;
    if fs3<fy
        fs3;
    else
        fs3=fy;
    end
    P(i)=(0.85*fc*a*b+As2*fs2+As3*fs3-As1*fs1)/1000
    M(i)=(0.85*fc*a*b*(h/2-a/2)+As3*fs3*(d-h/2)+As1*fs1*(d-h/2))/1000
    eb(i)=M(i)/P(i)      
    else
    disp('2.1) Acero en Tracción')
    As1=Asi; % Área en Tracción
    Es1=Eu*(d-C(i))/C(i);
    fs1=Es*Es1;1
    if fs1<fy
        fs1;
    else
        fs1=fy;
    end
    As2=Asm; % Área en Tracción
    Es2=Eu*(h/2-C(i))/C(i);
    fs2=Es*Es2;
    if fs2<fy
        fs2;
    else
        fs2=fy;
    end    
    disp('2.1) Acero en Compresión')
    As3=Ass; % Área en Compresión
    Es3=Eu*(C(i)-h+d)/C(i);
    fs3=Es*Es3;
    if fs3<fy
        fs3;
    else
        fs3=fy;
    end  
    P(i)=(0.85*fc*a*b+As3*fs3-As1*fs1-As2*fs2)/1000
    M(i)=(0.85*fc*a*b*(h/2-a/2)+As3*fs3*(d-h/2)+As1*fs1*(d-h/2))/1000
    eb(i)=M(i)/P(i)
    end
end
Pn=[P P0];
Mn=[M M0];
dint=plot(Mn,Pn,'marker','o'); % Grafica El diagrama de Interacción
text(Mb/3,1.5*Pb,'Falla a Compresión');
text(2*Mb/3,Pb/4,'Falla a Tracción');
xlabel('Mn(Ton-m)');
ylabel('Pn(Ton)');
title('DIAGRAMA DE INTERACCIÓN DE COLUMNA RECTANGULAR');
hold on
fbal=plot(Mb,Pb,'marker','o'); % Grafica el Punto de la Falla Balanceada
x=0:100:Mb;
y=(Pb/Mb)*x;
recta=plot(x,y); % Grafica la Recta de la Falla Balanceada
legend('Diagrama de Interacción','Falla Balanceada','Límite de Falla');
set(dint,'LineWidth',2,'LineStyle','-','Color','g','MarkerEdgeColor','r');
set(fbal,'LineWidth',2,'LineStyle',':','Color','b','MarkerFaceColor','b');
set(recta,'LineWidth',1,'LineStyle','-','Color','b');
grid on;
disp('FIN DE PROGRAMA')