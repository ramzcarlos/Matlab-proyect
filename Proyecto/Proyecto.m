clc
clear
% clearvars -except k k_v K M

%Maestria en Ciencias de la Ingeniría
%Sonia Valdez Mejia

% ----------- DATOS DE -----------------
% ------------ ENTRADA ------------------

% A partir de un archivo de texto
% load k_v.txt
% k=k_v;
% load K.txt
% load M.txt

% % A partir de valores numéricos
%solicita al usuario cuantos archivos va analizar
prompt = 'Numero de elementos del vector de rigidez: \n';
% lo guarda en la variable n_file
n_k = input(prompt);
k=zeros(1,n_k);
disp('******************************');
disp('');
disp('LLenado del vector de regidez');
disp(' ');
% Se define  el vector k
for i=1:n_k
    s1= "Ingresa el elemento ";
    s2= "  del vector de regidez\n";
    chr = int2str(i);
    s = strcat(s1,chr,s2);

    prompt = s;
    % lo guarda en la variable n_file
    k(i)=input(prompt);
end
disp('********************************** ');
%k=[30,70,50];
%k=k_v;

% Se define la matriz de rigidez K
disp('LLenado de la matriz de regidez ');

disp(' ');
K=zeros(n_k);
for i=1:n_k
    for j=1:n_k
    s1= "Ingresa el elemento ";
    s2= "  de la matriz  de regidez\n";
    chr = int2str(i);
    chr2 = int2str(j);
    s = strcat(s1,chr, chr2,s2);

    prompt = s;
    % lo guarda en la variable n_file
    K(i,j)=input(prompt);
    end
end
disp('****************************** ');

%K=[100,-70,0;-70,120,-50;0,-50,50];

% Se inicializa la matriz de masas M
disp(' LLenado de la matriz de masas');
disp(' ');
m_axu=zeros(1,n_k);
for i=1:n_k
    s1= "Ingresa el valor de  la masa ";
    s2= " (ton): \n";
    chr = int2str(i);
    s = strcat(s1,chr, s2);

    prompt = s;
    % lo guarda en la variable n_file
    m_aux(i)=input(prompt);
end

for i=1:n_k
    m_aux(i)=m_aux(i)/981;
end

M=diag(m_aux)

%
%M=zeros(n_k);
%for i=1:n_k
    %for j=1:n_k
    %s1= "Ingresa el elemento ";
    %s2= "  de la matriz  de masas\n";
    %chr = int2str(i);
    %chr2 = int2str(j);
    %s = strcat(s1,chr, chr2,s2);

    %prompt = s;
    % lo guarda en la variable n_file
    %M(i,j)=input(prompt);
    %end
%end
disp(' ');

%M=[0.1274,0,0;0,0.1019,0;0,0,0.2039];

% muestra en pantalla la información del vector k
fprintf('Vector de rigideces (ton/cm)\n');
disp(k);

% Muestra en pantalla la informacipon de la matriz K
fprintf('Matriz de rigidez (ton/cm)\n\n')
disp(K);

% Muestra en pantalla la información de la matriz M
fprintf('Matriz de masas (ton-s^2/cm)\n\n')
disp(M);


% Se guarda en la variable n el numero de elementos de la matriz K
n=length(K);

% Se crea un vector del 1, 2,3,..n elementos
Nivel=(1:n)';

% ****************************************************
% Inicio de los calculos de modos de vibración
%*****************************************************


% Cálculo de los valores y vectores característicos
% eig devuelve un vector de columna que contiene los valores propios de la
% matriz K y lo guarda en u
% eig devuelve un vector de columna que contiene los valores propios de la
% matriz M y lo guarda en d


[u,d]=eig(K,M);
fprintf('Eigen valores o valores característicos (w^2)\n\n');

%Genera un vector con los valores de la diagonal de la matriz d
Eigen_valores=diag(d);

% Frecuencias Modales
% disp('Frecuencias modales');
w=sqrt(Eigen_valores);

% Periodos Naturales
% disp('Periodos naturales');
T=(2*pi)./w;

disp('Periodos y frecuencias')
% Tabla_1=table(Nivel,w,T)
Modo_de_vibrar=(1:n)';
Tabla_1=table(Modo_de_vibrar,w,T);
disp(Tabla_1);

disp('w --> Frecuencias modales (rad/s)')
disp('T --> Periodos naturales (seg)')
disp(' ')

% Normalización de los vectores modales
% disp('Matriz modal normalizada')
modal=zeros(n);

for i=1:1:n
    for j=1:1:n
        modal(i,j)=u(i,j)/u(n,j);
    end
end
Matriz_modal_normalizada=table(Nivel,modal)

% ------------ANÁLISIS MODAL-----------------
% -------------------------------------------

% Registro de aceleraciones
% Vector de aceleraciones, para calculo de coeficientes de participación

%load sismo.txt

%*********************************************
disp('Lectura de archivos');
%ruta donde se encuentra el archivo txt y lo guarda en 
% 
    
prompt = 'Nombre del archivo para analizar: \n';
fileload= input(prompt,'s');
sismo=load(fileload);


z=input('Amortiguamiento modal (decimal) z= ');
disp(' ')

dt=input('Intervalo de tiempo en acelerograma (seg) dt= ');

Tabla_1=table(Modo_de_vibrar,T)

disp('El cálculo automático se basa en el método por carga lineal')
prompt=('Introducir manualmente las aceleraciones? Y/N\n');
Resp=input(prompt,"s");

if Resp=='Y' | Resp == 'y'
        % Crea un vector acc de 1, 2, 3, .. n posiciones y guardar ñps
        % valores de aceleración, modo acc
        acc=[1:n];
        for i = 1:n
        acc(i) = input(strcat('Ingrese la aceleración (cm/s^2) para T(',num2str(i),') del modo(',num2str(i),') acc(',num2str(i),')= '));
    end
    acc;
else

%%%%%%%%%%%%%%%%%%% FUNCIÓN POR CARGA LINEAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ingreso de aceleraciones tomadas del espectro
for i = 1:n
    acc(i) =Sa_T(sismo,z,T(i,1),1,dt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

disp(' ')
disp('Vector de aceleraciones por periodo (cm/s^2)')
Acel=acc';
Table_acc=table(Modo_de_vibrar,T,Acel)

% Vectores de modos de vibración
unos=ones(n,1);
for modo=1:n
    modoV{modo}=modal(:, modo);
end

% Coeficientes de participación por modo
cp=[1:n];

for i = 1:n
    cp(i) = ((transpose(modoV{i}))*M*unos)/(transpose(modoV{i})*M*modoV{i});
end

disp('Vector de factores de participación cp por modo')
disp(cp);
disp('Suma de coeficientes de participación')
SumaCP=sum(cp);
disp(SumaCP);
disp('OK! si la suma de los CPs es igual a 1.0')

% Máximas respuestas modales Desplazamientos Totales

for i=1:1:n
    for j=1:1:n;
        q_phi(i,j)=modal(i,j)*(cp(j)*acc(j)/w(j)^2);
    end
end

disp(' ')
disp('Matriz de desplazamientos totales modales por modo (cm)');
Matriz_q_phi=table(Nivel,q_phi);
disp(Matriz_q_phi);


% Desplazamientos relativos de entrepiso
%----------------------------------------

for i=2:1:n
    for j=1:1:n;
        q_Delta(1,j)=q_phi(1,j);
        for j=1:1:n;
            q_Delta(i,j)=q_phi(i,j)-q_phi(i-1,j);
        end
    end
end

disp('Desplazamientos relativos de entrepiso por modo (cm)')
Tabla_q_Delta=table(Nivel,q_Delta);
disp(Tabla_q_Delta);

% Cortante basal de entrepiso
%----------------------------

for i=1:1:n
    for j=1:1:n;
        Vb(i,j)=q_Delta(i,j)*k(i);
    end
end

disp('Matriz de cortantes de entrepiso por modo (Ton)')
Tabla_Vb=table(Nivel,Vb)

% Respuesta total SRSS
%------------------------

for i=1:1:n
    Modal_Vb(i,1)=sqrt(Vb(i,:)*Vb(i,:)');
end

% disp('Cortantes de entrepiso (SRSS)')
Modal_Vb;

for i=1:1:n
    Modal_Delta(i,1)=sqrt(q_Delta(i,:)*q_Delta(i,:)');
end

% disp('Desplazamiento relativo de entrepiso (SRSS)');
Modal_Delta;

for i=1:1:n
    Modal_u(i,1)=sqrt(q_phi(i,:)*q_phi(i,:)');
end

% disp('Desplazamiento total de entrepiso (SRSS)');
Modal_u;

% ------------ ANÁLISIS PASO A PASO -----------------
% Calculo de las ecuaciones modales; Metodo de paso a paso


% %%%%%%%%%%%%%%%%%%% FUNCIÓN POR CARGA LINEAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ingreso de aceleraciones tomadas del espectro
% Desde 1 a n llama la funcion Despazamiento,  Velocidad y Aceleracion
% manda como parametros el valor de cp*sismo, z, Tm dt
for i=1:1:n
    q_D(i,1)=DES(cp(i)*sismo,z,T(i),1,dt);
    q_V(i,1)=VEL(cp(i)*sismo,z,T(i),1,dt);
    q_A(i,1)=ACC(cp(i)*sismo,z,T(i),1,dt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cálculo de desplazamiento total (Paso a paso)
for i=1:1:n
    PaP_u(i,1)=(abs(modal(i,:)))*q_D;
end

%  disp('Desplazamiento total de entrepiso (Paso a paso)');
PaP_u;

% Cálculo de desplazamiento relativo (Paso a paso)
for i=1:1:n
    PaP_Delta(1)=PaP_u(1);
    for i=2:1:n;
        PaP_Delta(i,1)=PaP_u(i,1)-PaP_u(i-1,1);
    end
end

%  disp('Desplazamiento relativo de entrepiso (Paso a paso)');
PaP_Delta;

% Cálculo de cortante basal (Paso a paso)
for i=1:1:n
    PaP_Vb(i,1)=PaP_Delta(i)*k(i);
end

%  Cortantes de entrepiso (Paso a paso)
PaP_Vb;

disp('Máximos cortantes y desplazamientos modales (SRSS)')
Tabla_SRSS=table(Nivel,Modal_u,Modal_Delta,Modal_Vb)

disp('Modal_u     --> Desplazamiento total modal (cm)');
disp('Modal_Delta --> Desplazamiento relativo modal (cm)');
disp('Modal_Vb    --> Cortante modal de entrepiso y basal (Ton)');
disp(' ')

disp('Máximos cortantes y desplazamientos (Paso a paso)')
Tabla_Paso_a_Paso=table(Nivel,PaP_u,PaP_Delta,PaP_Vb)

disp('PaP_u     --> Desplazamiento total paso a paso (cm)')
disp('PaP_Delta --> Desplazamiento relativo paso a paso (cm)')
disp('PaP_Vb    --> Cortante paso a paso de entrepiso y basal (Ton)')
disp(' ')

for i=1:1:n
    Error_u(i,1)=(1-(Modal_u(i)/PaP_u(i)))*100;
    Error_Vb(i,1)=(1-(Modal_Vb(i)/PaP_Vb(i)))*100;
end

disp('Comparación y porcentaje de error')
Tabla_2=table(Nivel,Modal_u,PaP_u,Error_u,Modal_Delta,PaP_Delta,PaP_Vb,Modal_Vb,Error_Vb)

disp(' ')
disp('La variable "Error" está expresada en porcentaje % ')
disp(' ')
disp('               -----Variables-----')
disp('u     --> Desplazamiento total (cm)')
disp('Delta --> Desplazamiento relativo (cm)')
disp('Vb    --> Cortante de entrepiso y/o basal (Ton)')
disp(' ')
disp('               -----Prefijos-----')
disp('PaP   --> Calculado por método paso a paso')
disp('Modal --> Calculado por método modal espectral')

% --------- Gráfica de modos de vibrar ----------

for i=2:1:(n+1)
    for j=1:1:n
        Plot_modal(1,j)=0;
        Plot_modal(i,j)=modal((i-1),j);
    end
end
Plot_modal;
Resultados_modales=Plot_modal;
X=modal(:,1);
Y=[0:1:n];
Zero=zeros(n+1,1);
%Y es la aceleracion ya sea manual o automatica
for i=1:1:n
    subplot(1,n,i)
    plot(Resultados_modales(:,i),Y','*--',Zero,Y','*-k');
    title(['Modo' num2str(i)]);
    grid on
end


