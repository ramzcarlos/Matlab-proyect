%PROYECTO DINÁMICA ESTRUCTURAL, PASO A PASO
clc;
clear all;

%Introducir la matriz K, unidades en kg/m
k=[17595473.55	-16207217.74	4953526.626	-238807.7449	-185298.9556;-16207217.74	29322816.54	-18940129.61	3953704.988	-251689.4905;4953526.626	-18940129.61	27609544.38	-17294238.89	4008204.681;-238807.7449	3953704.988	-17294238.89	26323847.08	-12665728.61;-185298.9556	-251689.4905	4008204.681	-12665728.61	9051903.526]

%Introducir la matriz M, unidades kg s2/m
m=[2003.742359	0	0	0	0;0	1743.714724	0	0	0;0	0	1743.714724	0	0;0	0	0	1743.714724	0;0	0	0	0	1238.955199]


%Variables
C=[];
u=[];
phi=[];
suma=0;
vec=[];
Tp=[];

disp('Grados de libertad')
GDL=length(k)

%Vectores y valores característicos
Minv=inv(m);
A=Minv*k;
[u,d]=eig(A);

%Frecuencas modales
for i=1:1:GDL
    for j=1:1:GDL
        
        if(i==j)
             
            vec(i)=d(i,j);
            suma=suma+vec(i);
            W=vec';
        end
    end
end

disp('Frecuencia modal')
s=sqrt(W);
w=s(end:-1:1)

%Periodos
disp('Periodos modales 1-n')
for i=1:GDL
   T(i,:)=(2*pi)/w(i,1);
end
T

%Vector u
for i=1:GDL
   u(:,i)=(1/(u(1,i)))*u(:,i);
end
disp('Vector modal u_ij')
u;
u=fliplr(u)


%Decir respecto a que masa se debe normalizar
masa=input('Respecto a que masa desea normalizar ');
for i=1:GDL
   phi(:,i)=u(:,i)/u(masa,i);
end
disp('Normalizando vector u respecto a la masa n')
phi

NnAum=[0:1:GDL];
origAu=zeros(GDL+1,1);
orig=zeros(GDL,1);
for modo=1:GDL
    Modo{modo}=phi(:, modo);
end

MmodalAum=[orig';phi];
for modo=1:GDL
    X{modo}=MmodalAum(:, modo);
end

figure
          for kk=1:GDL
             subplot(1,GDL,kk)   
             plot(X{kk},NnAum','o--b',origAu,NnAum','o-k');
             ylabel('Nivel (n)');
             %ytickformat('%,f0')
             xlabel('Posición');
             title(['Modo' num2str(kk)]);
             sgtitle(['Formas modales']);    
             grid on 
          end

disp('Coeficientes de participación')
for i=1:GDL
    b=phi(:,i);
    VUnos=ones(GDL,1);
    cp_i(:,i)=(b'*m*VUnos)/(b'*m*b); %REVISAR EL NUMERO DE 1'S
end
cp_i

%***************************************************
    % Datos iniciales del método
    beta=1/4;
    gamma=1/2;
    desp_0=zeros(GDL,1); %u_0
    vel_0=zeros(GDL,1); %upunto_0
    phi_trans=transpose(phi);

    % solicita valor para el amortiguamiento en porcentaje y lo guarda en epsilon
    prompt = 'Amortiguamento (porcentaje) ';
    epsilon = input(prompt);
    % Lo convierte en porcentaje ejem 10 lo hace a 0.10
    epsilon = epsilon/100; 

%Cálculo de matriz de amortiguamiento
modo_13=(1/2)*[1/w(1,1) w(1,1);1/w(3,1) w(3,1)]
amort=[epsilon;epsilon];
a=inv(modo_13)*amort;
c=(a(1,1)*m)+(a(2,1)*k)

%solicita al usuario cuantos archivos va analizar
%prompt = 'Cuantos archivos deseas analizar: \n';
% lo guarda en la variable n_file
n_file = 1;
% crea un vector de tamaño de n_file ejem 1 o 2 o 3 dependiendo de lo que
% dio el usuario
if n_file == 1
    Name_sta=zeros(1:2);
else
    Name_sta=zeros(1:n_file);
end

Name_sta=Name_sta(:);
Name_sta=string(Name_sta);

    % Solicita valor de delta t, y lo almacena en at ejem 0.005
    prompt = 'valor de Delta t? ';
    at = input(prompt);

    %disp('Lectura de archivos');
    prompt = 'Nombre del archivo para analizar: \n';
    filename= input(prompt,'s');
    delimiterIn = ',';
    headerlinesIn = 1;
    As = load(filename);

    % Coloca los valores en una sola columna
    tam_as=size(As);
    as_i=tam_as(1);
    as_j=tam_as(2);
    tam_as=as_j*as_i;

    %Calculo de (qn)0 y (q.n)0
    for i=1:GDL
    q_n0(i,:)=(phi_trans(i,:)*m*desp_0)/(phi_trans(i,:)*m*phi(:,i));
    q_punto0(i,:)=(phi_trans(i,:)*m*vel_0)/(phi_trans(i,:)*m*phi(:,i));
    end

    %Cálculo de P0,M,C,K
    p0=-(m*VUnos*As(1,1));
    P0=phi_trans*p0;
    M=phi_trans*m*phi;
    C=phi_trans*c*phi;
    K=phi_trans*k*phi;

    %Cálculo de q..0
    acel_0=inv(M)*(P0-(C*q_punto0)-(K*q_n0));

    %Constantes a1,a2,a3
    a1=((1/(beta*at*at))*M)+(gamma/(beta*at))*C;
    a2=((1/(beta*at))*M)+((gamma/(beta))-15)*C;
    a3=(((1/(beta*2))-1)*M)+((gamma/(2*beta))-1)*C;

    %K con gorrito
    K_g=K+a1;
    
    %*********************************************
    %inicio de la segunda parte
    %*********************************************
    inv_K_g=inv(K_g);
    len_As=length(As);

    %empieza con el llanado de p
    p=zeros(len_As, GDL);
    
    for i=1:len_As-1
        p_aux=-m*VUnos*As(i+1);
        for j=1:GDL
            p(i+1,j)=p_aux(j);
        end
    end
    
    %comienza con el llenado de p_g
    %p_g (i+1) = (phi_trans * p (i+1)) + (a1*q_0 )+(a2*q_punto0)+(a3*acel_0)
    p_g=zeros(len_As, GDL);
    acel_0=acel_0(:);
    q=zeros(len_As, GDL);
    q_punto=zeros(len_As, GDL);
    q_dospuntos=zeros(len_As, GDL);
    for i=1:len_As
        if i>1
            p_aux=[];
            for j=1:GDL
                
                pg_a(j)=p(i,j);
                
                
            end
            
            
            pg_a=pg_a(:);
            
            

            %p_g_aux=(aux)+ (a1* q_n0)+(a2*q_punto0)+(a3*acel_0);
            p_g_aux=(phi_trans*pg_a)+ (a1* q(i-1))+(a2*q_punto(i-1))+(a3*acel_0);
            for j=1:GDL
                p_g(i,j)= p_g_aux(j);
                pq_a(j)=p_g(i,j);
            end
            pq_a=pq_a(:);
            %comienza el llenado de q
            %q (i+1) = la inversa de la matriz k_g * p_g(i+1)
            q_aux=inv_K_g*pq_a;
            for j=1:GDL
                q(i,j)= q_aux(j);
                qpunto_a(j)=q(i,j);
            end
            qpunto_a=qpunto_a(:);
            % comienza el llenado de q_punto
            %q_punto (i+1)   = ( gamma/(beta*at)) * (q(i+1)-q_n0) + (1-(gamma/beta)) * q_punto0 + at*(1-(gamma/(2*beta)) * acel_0
            
            q_punto_aux=(gamma/(beta*at)) * (qpunto_a-q(i-1)) + (1-(gamma/beta)) * q_punto(i-1) + at * (1-(gamma/(2*beta))) * acel_0;
            for j=1:GDL
                q_punto(i,j)= q_punto_aux(j);
                qdospuntos_a(j)=q(i,j);
            end
            qdospuntos_a=qdospuntos_a(:);
            %comeinza el llenado de q_dospuntos
            %q_dospuntos (i+1) = (1/(beta*(at^2)))*( q(i+1) – q_n0(i)) – (1/(beta*at))*q_punto0) – ((1/(2*beta))-1) * acel_0

            q_dospuntos_aux=(1/(beta*(at^2))) * (qdospuntos_a-q(i-1,j))-(1/(beta*at)*q_punto(i-1)) -((1/(2*beta))-1) * acel_0;

            
            for j=1:GDL
                
                q_dospuntos(i,j)=q_dospuntos_aux(j);
            end
        
        
        
        else
            for j=1: GDL
                
                q(i,j)=q_n0(j);
                p_g(i,j)=q_n0(j);
                q_punto(i,j)=q_punto0(j);
                q_dospuntos(i, j)= acel_0(j);
            end
    

        end
        
    end