clc;
clear all;

% valores constantes
fc=250;
fc1=fc*14.223343334285;
fy=4200;
b=30;
d=63;
e0=0.002;
r=5;
h=d+r;
Ec=14000*sqrt(fc);
Es=2080000;
As=6;
As1=As*2.87;
eu=0.00395631;
area=30.152916;
prompt= 'Se asume el valor de c de: \n';
c=input(prompt);
prompt= 'Se asume el valor de ec de: \n';
ec=input(prompt);
profundidad=[];

for i=1:11
    if i==1
        profundidad(i)=0;
    else
        profundidad(i)=profundidad(i-1)+(c/10);
    end
end

profundidad=profundidad(:);


brazo=[];
for i=1:11
    if i==11
        brazo(i)=0;
    else
        brazo(i)=((profundidad(i+1)-profundidad(i))/2)+profundidad(i);
    end
end

brazo=brazo(:);

deformacion=[];
for i=1:11
    deformacion(i)=((-profundidad(i)+c)*ec)/c;
end

deformacion=deformacion(:);

psi=[];
for i=1:11
    psi(i)=(((-0.5*fc1)/(eu-e0))*(deformacion(i)-eu))+(0.5*fc1);
end
psi=psi(:);

fc=[];
for i=1:11
    fc(i)=psi(i)*0.070307;
end
fc=fc(:);

F=[];
for i=1:11
    F(i)=area*fc(i);
end
F=F(:);

Fbrazo=[];
for i=1:11
    Fbrazo(i)=F(i)*brazo(i);

end
Fbrazo=Fbrazo(:);

Tabla=[profundidad brazo deformacion psi fc F Fbrazo];
FConcreto=sum(F);
FConcretob=sum(Fbrazo);
profundidadCc=FConcretob/FConcreto;

fprintf('\nValores de la tabla:\n\n');
fprintf('Profundad\tBrazo   \tDeformacion\tfc(psi)   \tfc (kg/cm2) \tF=A*fc, kg \tF*brazo\n');
for i=1:11
    for j=1:7
        fprintf('%f', Tabla(i,j));
        fprintf('\t');
    end
    fprintf('\n');
end

fprintf('Fuerza en el concreto= %f \t %f \n', FConcreto, FConcretob);
fprintf('Profumdidad de Cc= %f cm \n', profundidadCc);

