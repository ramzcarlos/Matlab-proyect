clear;
clc;
%***************************************************
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


for n=1: n_file
  
    % Solicita valor de delta t, y lo almacena en at ejem 0.005
    prompt = 'valor de Delta t? ';
    at = input(prompt);

    % solicita valor para el amortiguamiento en porcentaje y lo guarda en epsilon
    prompt = 'Amortiguamento (porcentaje) ';
    epsilon = input(prompt);
    % Lo convierte en porcentaje ejem 10 lo hace a 0.10
    epsilon = epsilon/100; 
    
    %*********************************************
    disp('Lectura de archivos');
    
    prompt = 'Nombre del archivo para analizar: \n';
    filename= input(prompt,'s');
    %Name_file(n)=filename;
    %filename = '04_dx.txt';
    delimiterIn = ',';
    headerlinesIn = 1;
    As = load(filename);
    % Coloca los valores en una sola columna
    tam_as=size(As);
    as_i=tam_as(1);
    as_j=tam_as(2);
    tam_as=as_j*as_i;

    %--------------------------
    %--------------------------
    %----------------------
    %solicita el valor inicial y lo guarda en t_i ejem 0.1
    %prompt = 'Periodo del oscilador: \n';
    T=[0.0745, 0.0918, 0.1258, 0.2013, 0.3617, 1.1646];
    tam_T=size(T);
    tam_T=tam_T(2);
    aux_T=tam_T+1;
    Resultados=zeros(aux_T, tam_as);
    for ini=1:tam_T
        %t_i = input(prompt);
        t_i = T(ini);
        t_f = t_i;
        t_c = t_i;
    
        % CRea el archivo de tiempo con los valores iniciales, incremento y finales
        Time=(t_i:t_c:t_f);
        % Colocal el vector Time en columna
        Time=Time(:);
        size_T=size(Time);
        size_T=size_T(1);

        Time=zeros(1,size_T+1);
        cont_incre=0;
        for i=1:size_T+1
            if i ==1
                Time(i)=0;
            elseif i ==2
                cont_incre=t_i;
                Time(i)=cont_incre;
            else
                cont_incre=cont_incre+t_c;
                Time(i)=cont_incre;
            end
        end

        Time=Time(:);

        % Asigna valores de beta y gamma
        beta=1/4;
        gamma=1/2;
    
        % s_t es el tamaño del vector T
        s_t=size(Time);
        %hold on
        % s_t lo convierte en entero 
        s_t=s_t(1);
        % se difine la matriz M que es donde se guardan todos los resultados
        M=zeros(s_t, n_file);
        % a_p es una variable de posicion  solo para el manejo de matriz
        a_p=1;
        %---------------------------



        %inicia el calculo por archivo es el ciclo que se va a repetir dependiendo
        %del numero de archivos a analizar

        %--------------------------
        %--------------------------


        aux_as=zeros(1,tam_as);
        cont_as=0;
        for i =1:length(As)
            for j =1:length(As(i,:))
                cont_as=cont_as+1;
                aux_as(cont_as)=As(i,j);
            end
        end
        As=aux_as;
        As = As(:);

        % Crear vector al tamaño de los archivos de entrada
        At=zeros(size(As));
        aux_at=0;
        %s gudara el tamaño del vector At
        s=size(At);
        % llena el vector a con inicio de 0 con incremento de delta t
        for i=1:s
            At(i)=aux_at;
            aux_at=aux_at+at;
        end
        % valor de la masa igual  a 1
        m=1;
    
        % Crea el vector Pi al tamaño del archivo de datos
        Pi=zeros(s);
    
        %llena el vector Pi con los datos de aceleracion de suelo multiplicado por
        %la masa
        for i=1:s
            Pi(i)=-m*As(i);
        end
    

    
        % Crea los vectores donde se colocaran los valores maximos por cada muestra
        % de tiempo t
        V=zeros(s_t);
        D=zeros(s_t);
        A=zeros(s_t);
        AT=zeros(s_t);

        %crea los vectores donde se iran almancenando los valores de cada corrida
        Acel=zeros(s);
        Psc=zeros(s);
        Desp=zeros(s);
        Vel=zeros(s);
        AcelTotal=zeros(s);
    
        % se empieza a mover en cada uno de los espacios del vector del 1 hasta el
        % tamaño del vector Time
        max_as=(max(abs(As)));
        for j=1:s_t
            if j ==1
                A(j)=(max_as);
                D(j)=(0);
                V(j)=(0);
                AT(j)=(max_as);
            else
                % Se asigan el primer valor del vector t que es el valor inicial que
                % indico el usuario
                t=Time(j);
        
                k=((2*pi)/t)^2;
                c=2*epsilon*sqrt(k*m);
                a1=((1/(beta*at*at))*m)+((gamma/(beta*at))*c);
                a2=((1/(beta*at))*m)+((gamma/(beta)-1)*c);
                a3=((1/(2*beta)-1)*(m))+(at*(gamma/(2*beta)-1)*c);
                omega=sqrt((k/m));
      
                % k1 es k con gorrito
                k1=k+a1;
    
                % Comienza a realizar los calculos para vector velocidad, desplazamiento,
                % aceleración y aceleración total
                for i=1:s
        
                    if i ==1
                        %se inicializan los ventores con los datos de inicio
                        Acel(i)=Pi(i);
                        Psc(i)=Pi(i);
                        Desp(i)=0;
                        Vel(i)=0;
                        AcelTotal(i)=Acel(i)+As(i);
                    else
                        %realiza los calculos y va llenando cada uno de los vectores
                        Psc(i)=Pi(i)+a1*Desp(i-1)+a2*Vel(i-1)+a3*Acel(i-1);
                        Desp(i)=Psc(i)/k1;
                        Vel(i)=(2/at)*(Desp(i)-Desp(i-1))-Vel(i-1);
                        Acel(i)=((4/(at*at))*(Desp(i)-Desp(i-1)))-((4/at)*Vel(i-1))-Acel(i-1);
                        AcelTotal(i)=Acel(i)+As(i);
                    end
                end
                % Busca el maximo de los los vectores en valor absouto y lo almancena
                % en el espacio que corresponde al valor Time que se este realizando la
                % corrida
                A(j)=(max(abs(Acel)));
                D(j)=(max(abs(Desp)));
                V(j)=(max(abs(Vel)));
                AT(j)=(max(abs(AcelTotal)));
            end
        end
       
        %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         % Asigna los 4 vectores(desplazamiento, velocidad, aceleracion y
         % aceleracion total) en la matriz M
        for zi=1:length(D)
            M(zi, a_p)=D(zi);
        end
        a_p=a_p+1;
        for zi=1:length(V)
            M(zi, a_p)=V(zi);
        end
        a_p=a_p+1;
        for zi=1:length(A)
            M(zi, a_p)=A(zi);
        end
        a_p=a_p+1;
        for zi=1:length(AT)
            M(zi, a_p)=AT(zi);
        end
        a_p=a_p+1;

             
            
    for fi=1:tam_as
        Resultados(ini, fi)=Desp(fi);
    end
        



    
    end
end

for f=1:tam_as
    suma_r=0;
    for fi=1:ini
        suma_r=suma_r+Resultados(fi, f);
    end
    Resultados(aux_T, f)=suma_r;
end

for i=1:aux_T
    graf=zeros(1, tam_as);
    for j=1:tam_as
        graf(j)=Resultados(i,j);
    end
    c_plot=string(i);
    f=strcat("f",c_plot);
    f=figure
    if i==7
    tte=strcat('Desplazamiento relativo:   ', 'sumatoria');
    else
      tte=strcat('Desplazamiento relativo:   ', c_plot);  
    end
    plot(At,graf,'k')
         %el titulo que tenga el tte lo pone en la grafica
         title(tte)
         xlabel("T (s)")
         ylabel("Desplazamiento (cm)")
end
