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




%----------------------
%solicita el valor inicial y lo guarda en t_i ejem 0.1
prompt = 'Periodo del oscilador: \n';
t_i = input(prompt);
t_f = t_i;
t_c = t_i;
    
% CRea el archivo de tiempo con los valores iniciales, incremento y finales
T=(t_i:t_c:t_f);
% Colocal el vector T en columna
T=T(:);
size_T=size(T);
size_T=size_T(1);

T=zeros(1,size_T+1);
cont_incre=0;
for i=1:size_T+1
    if i ==1
        T(i)=0;
    elseif i ==2
        cont_incre=t_i;
        T(i)=cont_incre;
    else
        cont_incre=cont_incre+t_c;
        T(i)=cont_incre;
    end
end

T=T(:);

% Asigna valores de beta y gamma
beta=1/4;
gamma=1/2;
    
% s_t es el tamaño del vector T
s_t=size(T);
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
for n=1: n_file
    % solicita al usuario datos de entrada
    %prompt = 'Nombre de la estacion: \n';
    %sta = input(prompt,'s');
    
    %prompt = 'Pais: \n';
    %contry = input(prompt,'s');

    %prompt = 'Ciudad o Estado: \n';
    %town = input(prompt,'s');
    %prompt = 'fecha: \n';
    %date = input(prompt,'s');
    %info=strcat(contry,"-", town);
   
    %Name_sta(n)=info;
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
    %ruta donde se encuentra el archivo txt y lo guarda en 
    % 
    
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
    
    %Energias acumuladas

    E_kin_acum = zeros(s);
    E_strain_acum = zeros(s);
    E_damp_acum = zeros(s);
    E_total_acum = zeros(s);
    E_kin_sa = zeros(s);
    E_strain_sa = zeros(s);
    E_damp_sa = zeros(s);
    E_total_sa = zeros(s);
    E_input_acum = zeros(s);
    E_input_sa = zeros(s);

    %crea los vectores donde se iran almancenando los valores de cada corrida
    Acel=zeros(s);
    Psc=zeros(s);
    Desp=zeros(s);
    Vel=zeros(s);
    AcelTotal=zeros(s);
    
    % se empieza a mover en cada uno de los espacios del vector del 1 hasta el
    % tamaño del vector T
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
            t=T(j);
        
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
            % en el espacio que corresponde al valor T que se este realizando la
            % corrida
            A(j)=(max(abs(Acel)));
            D(j)=(max(abs(Desp)));
            V(j)=(max(abs(Vel)));
            AT(j)=(max(abs(AcelTotal)));
        end
    end
    
    %termina el programa de hacer todos los calculos
    %*************************************************************************
    % Seccion para guardar los resultados en un archivo txt
    % solicita al usuario el nombre con el que se guardan los resultados
    %prompt = 'Guardar los resultados de manera individual? Y/N [Y]: ';
    %resp = input(prompt,'s');
    %%max_as=(max(abs(As)));
    %if resp == 'Y'| resp == 'y'
       
    


        %prompt = 'Nombre del archivo para guardar resultados para el trabajo \n';
        %disp(n)
        %file_salve = input(prompt,'s');
        %file_path_salve=pwd+"/"+file_salve;
        %fileID = fopen(file_path_salve,'w');
        %formatSpec = '%.3f\t%.3f\t%.3f\t%.3f';
        %fprintf(fileID,"nombre de la estacion:  "+sta+"\n");
        %fprintf(fileID,"Pais: "+contry+"\n");
        %fprintf(fileID,"Ciudad o Estado: "+town+"\n");
        %fprintf(fileID,"Fecha: "+date+"\n");
        %fprintf(fileID,"Delta t: "+at+"\n");
        %fprintf(fileID,"Valor maximo de aceleracion de suelo: "+max_as+"\n");
        %fprintf(fileID,"id\tDes\tVel.\tAcel. Total.\n");
        %for data =1:length(A)
             %%fprintf("%d %d \n", i, j)
            %fprintf(fileID, formatSpec,T(data), D(data), V(data), AT(data));
            %fprintf(fileID,"\n");
       % end
      %% fprintf(fileID,"\n");
    
        %fclose(fileID);
    %end
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

             
            

        



    

end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%prompt = 'Guardar los resultados de manera grupal? Y/N [Y]: ';
%resp = input(prompt,'s');
%if resp == 'Y'| resp == 'y'
     %prompt = 'Nombre del archivo para guardar resultados para el trabajo \n';
     %disp(n)
     %file_salve = input(prompt,'s');
     %file_path_salve=pwd+"/"+file_salve;
     %fileID = fopen(file_path_salve,'w');
     %formatSpec = '%f %f';
     %tam=size(M);
     %n=tam(1);
     %m=tam(2);
    %for i =1:n
        %for j =1:m
            %fprintf(fileID, formatSpec, M(i,j));
         %end
        %fprintf(fileID,"\n");
    %end
     % fprintf(fileID,"\n");
    %fclose(fileID);
%end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


cont_plot=1;
fprintf("\n");
%prompt = 'Deseas realizar una grafica de manera individual? Y/N [Y]: ';
%resp = input(prompt,'s');
resp='y'
while resp == 'Y'  | resp == 'y'
    if n_file ==1
        datos=1;
    else

        %prompt = '\n Indique el de datos quieres graficar?  ';
        %datos = input(prompt);
        datos=1
        while datos < 0 & datos > n_file
            disp("Datos fuera de rango");
            prompt = '\nIndique el conjunto de datos que quieres graficar?  ';
            datos = input(prompt);
        end
    end
    fprintf("\n")
    %prompt = 'Grafica desplazamiento (D), velocidad (V), Aceleracion Total (A) ';
    %grafica = input(prompt, "s");

    grafica='D'

    while grafica ~= "D" & grafica ~= "V"& grafica ~= "A"
        disp ("opción incorrecta por favor indica que tipo de gafica deas desarrollar");
        prompt = 'Grafica desplazamiento (D), velocidad (V), Aceleracion Total (A) ';
        grafica = input(prompt, "s");
    end
    n=1;
    name_plot=Name_sta(datos);
    if n_file ==1
        g_d=1;
        g_v=2;
        g_a=4;
    else
        for i=2:datos
            n=n+4;
        end
        g_d=n;
        g_v=n+1;
        g_a=n+3;
    end
    %hold on
    if grafica =="D"
         c_plot=string(cont_plot);
         f=strcat("f",c_plot);
         f=figure
         %gráfica del sismo

         %aqui es donde se pone el titulo de las graficas y se cuarda en en
         %la variable tte
         tte=strcat('Desplazamiento relativo:   ', name_plot);
         plot(At,Desp,'k')
         %el titulo que tenga el tte lo pone en la grafica
         title(tte)
         xlabel("T (s)")
         ylabel("Desplazamiento (cm)")
         
         cont_plot=cont_plot+1;
    elseif grafica =="V"
         c_plot=string(cont_plot);
         f=strcat("f",c_plot);
         f=figure
        
       
         
         %xlim([0 4])plot(T,M(1:s_t,g_v))
         tte=strcat('Velocidad relativa:   ', name_plot);
         plot(At,Vel,'k')
         title(tte)
         xlabel("T (s)")
         ylabel("Velocidad (cm/s)")

         %hold off
         cont_plot=cont_plot+1;

    else
        %hold on
         c_plot=string(cont_plot);
         f=strcat("f",c_plot);
         f=figure
  
        
        tte=strcat('Aceleracion total:   ', name_plot);
        plot(At,AcelTotal,'k')
         title(tte)
         xlabel("T (s)")
         ylabel("Aceleración (cm/s^2)")

         cont_plot=cont_plot+1;
        %hold off
    end
    
    %fprintf("***************************************\n");
    %fprintf("Usted a creado la grafica: %s\n", tte);
    %fprintf("***************************************\n\n")

    %prompt = 'Deseas realizar otra grafica de manera individual? Y/N [Y]: ';
    %resp = input(prompt,'s');
    resp = 'n';
    fprintf("\n");
end
 
%---------------------------------------------------------------------

cont_plot=1;
if n_file > 1
    prompt = 'Deseas realizar una grafica de manera grupal? Y/N [Y]: ';
    resp = input(prompt,'s');
    fprintf("\n");
	while resp == 'Y'  | resp == 'y'
    
    	prompt = 'Grafica desplazamiento (D), velocidad (V), Aceleracion (A) ';
    	grafica = input(prompt, "s");

    	while grafica ~= "D" & grafica ~= "V"& grafica ~= "A"
        	disp ("opción incorrecta por favor indica que tipo de gafica deas desarrollar");
        	prompt = 'Grafica desplazamiento (D), velocidad (v), Aceleracion (A) ';
        	grafica = input(prompt, "s");
    	end
    	n=1;

    	dis_grafica=["--k", ":k","-.k","-k", "-xk", ":ok"];

    	
        g_d=1;
        g_v=2;
        g_a=4;

    	if grafica =="D"
         
         	f1=figure
         	hold on
         	for i=1:n_file
             	plot(T,M(1:s_t,g_d), dis_grafica(i),'Linewidth',0.5 )
             	ttg="Desplazamiento relativo"
             	title(ttg)
             	xlabel("T (s)")
             	ylabel("Desplazamiento (cm)")
             	g_d=g_d+4;
             	legend(Name_sta)
         	end
         	hold off
         
    	elseif grafica =="V"
        	f2=figure
        	hold on
        	for i=1:n_file
            	plot(T,M(1:s_t,g_v), dis_grafica(i), 'Linewidth',0.5)
            	ttg="Velocidad relativa";
            	title(ttg)
            	xlabel("T (s)")
            	ylabel("Velocidad (cm/s)")
            	g_v=g_v+4;
            	legend(Name_sta)
        	end

        	hold off
         	

    	else
        	
       
        	f3=figure

        	
        	hold on
        	for i=1:n_file
            	plot(T,M(1:s_t,g_a), dis_grafica(i),'Linewidth',0.5 )
            	ttg="Aceleracion Total";
            	title(ttg)
            	xlabel("T (s)")
            	ylabel("Aceleración (cm/s^2)")
            	g_a=g_a+4;
            	legend(Name_sta)
        	end
        	
        
        	hold off
    	end
    

    
        fprintf("***************************************\n");
        fprintf("Usted a creado la grafica grupal de: %s\n", ttg);
        fprintf("***************************************\n\n")
    	prompt = 'Deseas realizar otra grafica de forma grupal? Y/N [Y]: ';
    	resp = input(prompt,'s');
        fprintf("\n");
	end
end
%gráfica  del acelerograma
figure (4)
tte=strcat('Acelerograma:  ');
         plot(At,As,'k')
         %el titulo que tenga el tte lo pone en la grafica
         title(tte)
         xlabel("T (s)")
         ylabel("Aceleración (cm/s^2)")
% CALCULO DE ENERGÍAS

for j=2:tam_as
    E_kin = m*Vel(j)*Vel(j)/2.0;
    E_kin_acum(j) = E_kin_acum(j-1) + E_kin;
    E_kin_sa(j) = E_kin;
    E_strain = (1/2)*k*Desp(j)*Desp(j);
    E_strain_acum(j) = E_strain_acum(j-1) + E_strain;
    E_strain_sa(j) = E_strain;
    E_damp = c*Vel(j)*at;
    E_damp_acum(j) = E_damp_acum(j-1) + E_damp;
    E_damp_sa(j) = E_damp;
    E_total_acum(j) = E_kin_acum(j) + E_strain_acum(j) + E_damp_acum(j);
    E_total_sa(j) = E_kin_sa(j) + E_strain_sa(j) + E_damp_sa(j);
    E_input = -m*As(j)*(Desp(j)-Desp(j-1));
    E_input_acum(j) = E_input_acum(j-1) + E_input;
    E_input_sa(j) = E_input;

end

        figure (5)
        %plot(At,E_kin_acum,"--k")
       
        plot(At,E_kin_acum,"--", At, E_strain_acum, ":", At,E_total_acum,"-", At,E_damp_acum,"-." )
         title('Energía disipada');
        xlabel("Tiempo (s)");
        ylabel("Energía (cm/s)^2");
        %colororder(["#8040E6";"#1AA640";"#E68000"])

        legend({'Energía cinética','Energía de deformación', 'Energía total', 'Energía de amortiguamiento'},'Location','northwest')
        %hold on

        figure (6)
        
        plot(At,E_kin_sa,"--", At, E_strain_sa, ":", At,E_total_sa,"-", At,E_damp_sa,"-.",At,E_input_acum,"-x")
        title('Energía disipada sin acumular');
        xlabel("Tiempo (s)");
        ylabel("Energía (cm/s)^2");
        %colororder(["#8040E6";"#1AA640";"#E68000"])

        legend({'Energía cinética','Energía de deformación', 'Energía total', 'Energía de amortiguamiento','Energía de entrada'},'Location','northwest')
         %   plot(At,E_strain_acum,":k")
         %   legend('Energía de deformación')
         %   plot(At,E_total_acum,"-k")
         %   legend('Energía total')
         %   plot(At,E_damp_acum,"-.k")
         %   legend('Energía de amortiguamiento')

        %hold off



%genera las graficas en una matriz de 2 * 2

%{   

t = tiledlayout(2,2);
ax1=nexttile
hold on
axi=1;
for i=1:n_file    
    hold on
    plot(ax1,T,M(1:s_t,axi))
    title(ax1,'desplazamiento relativo')
    axi=axi+4
    %xlabel(t,'x-values')
    %ylabel(t,'y-values')
end   
hold off
ax2=nexttile
hold on
axi=2;
for i=1:n_file
    %hold on

    plot(ax2,T,M(1:s_t,axi))
    axi=axi+4
    title(ax2,'velocidad relativa')
end
%hold off

    
ax3=nexttile
%xlabel(t,'x-values')
%ylabel(t,'y-values')
hold on
axi=3;
for i=1:n_file 
    hold on
    plot(ax3,T,M(1:s_t, axi))
    axi=axi+4
    title(ax3, 'aceleración relativa')
    %xlabel(t,'x-values')
    %ylabel(t,'y-values')
end
hold off
ax4=nexttile
hold on
axi=4;
for i =1: n_file
    %hold on
    plot(ax4,T,M(1:s_t, axi))
    axi=axi+4;
    title(ax4,'aceleración total')
    %xlabel(t,'x-values')
    %ylabel(t,'y-values')
end
hold off
%}
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




