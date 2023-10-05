%% Parámetros del código

% Definimos las variables del material de los engranajes globales
global HB St rho_eng S_ut rho_eje a Vel_max W Rpm_max Rpm_min Tor_max Potencia rueda_carga R_rodadura R_pendiente R_inercia R_aerodinamica

%% Constantes del material engranajes  " AISI 8620 " http://ferrumaceros.com/assets/files/Ficha-AISI8620.pdf
rho_eng=7.87*1000; %Densidad del material [g/cm^3] a [kg/m^3]
%k_adm=80         %Acero aleado, cementado y templado [kg/cm^2]
HB=670;           %Dureza Brinnel

%Resistencia a flexión [Psi] a [MPa]
%St=(-274+167*HB-0.152*HB^2)/145;  %Máximo grado 1
St=(6235+174*HB-0.126*HB^2)/145;   %Máximo grado 2

%% Constantes del material eje " AISI 4140 " https://www.cga.com.co/wp-content/uploads/2020/07/Ficha_T%C3%A9cnica_Aceros_Grado_Ingenier%C3%ADa_4140.pdf
sigma_y=75*9.807;   %Limite fluencia [kgf/mm^2] a [MPa]
S_ut=88*9.807;      %Resistencia traccion [kgf/mm^2] a [MPa]
rho_eje=7.87*1000;  %Densidad del material [g/cm^3] a [kg/m^3]

%% Parametros del vehiculo
% Definimos las variables del vehiculo globales
a=50/(3.6*8);  %Aceleración: 0 a 50 [km/h] en 8 segundos
Vel_max=80;    %Velocidad máxima [km/h]
W = 750*9.81;  %Peso total del vehiculo en [N]

% Llanta 125/80 R13
f=0.9; %Factor de compensación llanta bajo carga
%rueda_carga=f*(13*25.4 +0.8*125*2)/(2*1000); %Radio de la llanta bajo carga
rueda_carga=0.265;

%% Parametros de los motores
% Cargar tabla información motores
% Estructura del archivo Infomotores
% || Id motor || Rpm max || Rpm min || Par max || Potencia ||
load('Motores.mat');

Num_motor=2; %Fila de datos del motor a usar según la tabla InfoMotores.mat

% Definimos las variables del motor globales
Tor_max=Motores(Num_motor,4);   %Torque max motor [Nm]
Potencia=Motores(Num_motor,5);  %Potencia maxima motor [kW]
Rpm_max=Motores(Num_motor,2);   %Rpm máximas del motor 
Rpm_min=Motores(Num_motor,3);   %Rpm mínimas del motor

%% Cálculos - llamado de funciones %%

% Calculo de resistencias y potencia mínima
[R_rodadura, R_pendiente, R_aerodinamica, R_inercia, P_min]=Calculo_Pot_min(Vel_max, W, a);

% Estructura de la matriz de los intervalos de las relaciónes
% || ID || Primera min || Primera_max || Segunda min || Segunda max
Intervalos=Calculo_Relaciones(Motores, Vel_max, rueda_carga, W, R_rodadura);

% Optimización
fun=@moo_functions; % Llamada del archivo con las funciones objetivo
nonlcon=@moo_const; % Llamada del archivo con las restricciones

nvars=2; % Número de variables
A=[]; b=[]; Aeq=[]; beq=[];

% Definimos los intervalos de las variables
% Cada limite se define el mismo orden en que se asignan las variables en
% los otros archivos a llamar
lb=[Intervalos(Num_motor,2), Intervalos(Num_motor,4)]; %Se definen los limites inferiores
ub=[Intervalos(Num_motor,3), Intervalos(Num_motor,5)]; %Se definen los limites superiores

[x, fval, exitflag, output] = gamultiobj(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon);

% Al momento de graficar, multiplicamos los resultados de las 
% funciones a maximizar por menos para obtener el valor real
Peso=fval(:,1);    %Peso
Torque=-fval(:,2); %Torque

I_primera=x(:,1);
I_segunda=x(:,2);

T=table(Peso,Torque, I_primera, I_segunda); %Tabla con las variables de interés
T_sorted=sortrows(T,{'I_primera'})       %Tabla excluyendo valores repetidos de la primera relación
%writetable(T_sorted, 'Optimizacion.xlsx');

Peso_sorted=T_sorted.Peso;     %Se obtienen los datos del peso
Torque_sorted=T_sorted.Torque; %Se obtienen los datos del torque

plot(Peso_sorted, Torque_sorted, '*-')
xlabel('Peso [N]')
ylabel('Torque [Nm]')
grid on

%% Resultados (solución escogida)

I_primera_f=1.5981; I_segunda_f=0.74252;
[z11, z12, z21, z22, m, U1, Fr1, Fa1, U2, Fr2, Fa2]=Calculo_engranajes(I_primera_f, I_segunda_f, HB, Rpm_max, Tor_max, St);

[E1]=Esp_engranajes(z11, m); [E2]=Esp_engranajes(z12, m); [E3]=Esp_engranajes(z21, m); [E4]=Esp_engranajes(z22, m);
Datos_engranajes=[E1', E2', E3', E4']
%writematrix(Datos_engranajes, 'Datos engranajes.xlsx');

[d1e, R1_Ay_1e, R1_Az_1e, R2_Ay_1e, R2_Az_1e, R1_Bx_1e, R1_By_1e, R1_Bz_1e, R2_Bx_1e, R2_By_1e, R2_Bz_1e]=Eje_principal(U1, Fr1, Fa1, z11, U2, Fr2, Fa2, z21, m, S_ut);
[d2e, R1_Ax_2e, R1_Ay_2e, R1_Az_2e, R2_Ax_2e, R2_Ay_2e, R2_Az_2e, R1_By_2e, R1_Bz_2e, R2_By_2e, R2_Bz_2e]=Eje_secundario(U1, Fr1, Fa1, z12, U2, Fr2, Fa2, z22, m, S_ut);

rodamientos=[calc_rod_bolas(R1_Ay_1e, R1_Az_1e, R2_Ay_1e, R2_Az_1e, Rpm_max), calc_rod_rodillos(R1_Bx_1e, R1_By_1e, R1_Bz_1e, R2_Bx_1e, R2_By_1e, R2_Bz_1e, Rpm_max);
    calc_rod_bolas(R1_By_2e, R1_Bz_2e, R2_By_2e, R2_Bz_2e, Rpm_max), calc_rod_rodillos(R1_Ax_2e, R1_Ay_2e, R1_Az_2e, R2_Ax_2e, R2_Ay_2e, R2_Az_2e, Rpm_max)];

% Rodamientos seleccionados:
% Bolas 1 eje: SKF 62207-2RS1
% Agujas 1 eje: SKF K 35x45x20
% Bolas 2 eje: SKF 6305
% Agujas 2 eje: SKF K 25x32x16

%Componentes=Calculo_embrague(Motores, Num_motor);

%% Planteamiento de funciones a utilizar %%

% Resistencias y potencia mínima [KW]
function  [Fr, Fp, Fa, Fj, Pot_min]=Calculo_Pot_min(Vel_max, Peso, Aceleracion)
    
    %Las resistencias se requieren para calcular la potencia mínima
    %requerida en el vehiculo y el torque generado en las ruedas

    V = Vel_max/3.6; %Velocidad en [km/h] a [m/s]
    
    %% Resistencia a la rodadura
    %La velocidad debe estar en [mph]
    fo = 0.009;  
    fs = 0.004;
    Cfr = fo +  3.24*fs*((V*2.237)/100)^2.5; %coeficiente de resistencia a la rodadura 
    Fr = Cfr*Peso;
    
    %% Resistencia aerodinamica
    CD = 0.39;      %coeficiente de resistencia aerodinamica
    A = 1.40875;    %area frontal del vehiculo [m^2] 
    rho_a = 1.204;  %densidad del aire a T ambiente [kg/m^3]
    Fa = CD*A*(rho_a/2)*V^2;

    %% Resistencia pendiente
    pendiente=0.26; %Pendiente avance vertical por cada 100 metros horizontales X/100 (16°)
    Fp=Peso*pendiente*0.375;

    %% Resistencia inercia
    j=Aceleracion;
    Fj=(Peso/9.81)*j;
    
    %% Potencia mínima
    Pot_min = ((Fa+Fr+Fj)*V/1000)*0.6;
    % Se multiplica por 0.6 debido a que la velocidad usada en las 3 fuerzas fue de 80 km/h
    % pero en los casos analizados se asume que el vehiculo va a 50 km/h
end

% Intervalos relaciones de la caja
function cambios=Calculo_Relaciones(Motores, Vel_max, R_rueda_carga, Peso, Fr)

    Dimension_matriz=size(Motores);  %Motores es el archivo "InfoMotores"
    cambios=zeros(Dimension_matriz(1),5);
    
    %% Constantes
    N_eficiencia=0.96; %Factor seguridad por perdidas de eficiencia
    miu_coeff = 0.8;   %Coeficiente de friccion estatica en condiciones de superficie seca a 90m/h
    alpha_st = 15.1;   %Angulo máximo de inclinación de la carretera
    
    I_E=3.7;          %Constante relación transmisión final-optimizar torque
    
    for i= 1:Dimension_matriz(1)
        
        %Ecuaciones de: 
        % "Design and Optimization of a Drivetrain with Two-speed Transmission for Electric Delivery Step Van"
        I_primera_min= (R_rueda_carga*Peso*((Fr/Peso)*cosd(alpha_st)+sind(alpha_st)))/(Motores(i,4)*N_eficiencia*I_E);
                        
        I_segunda_max= (3.6*pi*R_rueda_carga*Motores(i,2))/(30*Vel_max*I_E);
        I_segunda_min= (3.6*pi*R_rueda_carga*Motores(i,3))/(30*Vel_max*I_E);
        
        %Ecuacion de:
        % "Comprehensive design and optimization of an electric vehicle powertrain equipped with a two-speed dual-clutch transmission"
        I_primera_max=(R_rueda_carga*Peso*miu_coeff)/(Motores(i,4)*N_eficiencia); %Se debe dividir por I_E?
        
        cambios(i,1)=Motores(i,1); 
        cambios(i,2)=I_primera_min; cambios(i,4)=I_segunda_min;
        cambios(i,3)=I_primera_max; cambios(i,5)=I_segunda_max;
    end
    
end

% Calculo de fuerzas,dientes y modulo de los engranajes
function [z11, z12, z21, z22, m, sigma_f1, sigma_f2, U1, Fr1, Fa1, U2, Fr2, Fa2]=Calculo_engranajes( ...
    I_primera, I_segunda, HB, Rpm_max, Tor_max, St)

% Número de dientes y relaciones reales
% Consideraciones:
% Todos engranajes misma distancia entre centros
% Engranajes con el mismo módulo
% Constante la suma de número de dientes de cada par de engranajes
% Las relaciones se ingresan en valores >1 
    
    beta=23;    %Angulo de helice
    alpha=20;   %Angulo de presión
    %zn=14;      %Zn >=14 número de dientes virtual 

    NS=1.5;

    % Primer conjunto
    %z11=ceil(zn*cosd(beta)^3); %Número de dientes real - se escoge el mínimo
    z11=16;
    z12=ceil(z11*I_primera);
    zt=z11+z12;
    i1=z12/z11;                %Primera relación dada por los números de dientes
    
    % Segundo conjunto
    %z21=round(cons/(1+(1/I_segunda)));
    z21=round(zt/(1+I_segunda));
    z22=zt-z21;
    i2=z22/z21; %Segunda relación dada por los números de dientes

    % Módulo y distancia entre ejes
    
    T_vida=250000; %Tiempo de vida caja estimado [km]
    Vel_media=60;  %Velocidad media del vehiculo [km/h]
    E=2100000;     %Módulo elasticidad acero [kg/cm^2] 
    
    T_duracion=T_vida/Vel_media; %Horas duración caja

    t=(Rpm_max*T_duracion*60)/10^6;
    k_adm=6800*HB^2/(t^(1/3)*E);
        
    load('VidaUtil.mat') %Se carga el archivo de vida util engranajes segun porcentaje
    Dimension_matriz=size(Vida_util);

    for i= 1:Dimension_matriz(1)
        Vida_util(i,3)=T_duracion*Vida_util(i,2); %Horas utiles de cada marcha
        phi=16.895*Vida_util(i,3)^(-0.332);
        Vida_util(i,4)=phi*k_adm;
    end

    psi=12;              %Factor guiado asumiendo calidad y condiciones normales
    Tor=Tor_max*10.1972; %Conversion KW a CV

    % Primer conjunto - modulo en [mm]
    m1=(2*Tor*(i1+1)*(cosd(beta)^4)/(Vida_util(1,4)*z11^2*psi*cosd(alpha)*sind(alpha)*i1))^(1/3)*10; %Fallo superficial

    % Segundo conjunto - modulo en [mm]
    % El modulo se halla a partir de los dientes del piñon
    if (z21<z22)
        zmin=z21; % El modulo se halla a partir de los dientes del piñon
    else
        zmin=z22;
    end
    
    m2=(2*Tor*(i2+1)*(cosd(beta)^4)/(Vida_util(2,4)*zmin^2*psi*cosd(alpha)*sind(alpha)*i2))^(1/3)*10; %Fallo superficial
    
    if (m1>m2) % Se escoge el modulo mayor - se aproxima entero mayor
        m_round=round(m1);

        if (m1+0.5==m_round)
            m=m1;
        elseif (m_round>m1)
            m=m_round;
        else 
            m=m_round+0.5;
        end
    else
        m_round=round(m2);

        if (m2+0.5==m_round)
            m=m2;
        elseif (m_round>m2)
            m=m_round;
        else
            m=m_round+0.5;
        end
    end
    
    sigma_adm=0;
    sigma_f2=1;
    sigma_f1=1;

    while (sigma_adm/sigma_f1<NS)

    b=m*psi; % Ancho de diente
    
    %% Fuerzas y comprobación sobre los engranajes 
    
    T_max=Tor_max*1000; %Torque maximo [Nmm]
    
    % Primera marcha
    R_p11=(m/2)*(z11/cosd(beta)); %Radio primitivo piñon
    U1=T_max/R_p11;               %Fuerza tangencial
    Fr1=U1*tand(alpha);
    Fa1=U1*tand(beta);
    %W1=U1/(cosd(beta)*cosd(alpha));

    % Segunda marcha
    R_p21=(m/2)*(zmin/cosd(beta)); %Radio primitivo piñon  
    U2=T_max/R_p21;                %Fuerza tangencial
    Fr2=U2*tand(alpha);
    Fa2=U2*tand(beta);
    %W2=U2/(cosd(beta)*cosd(alpha));

    %% Factores comportamiento

    % Factor aplicacion, fuente de energia uniforme y carga con impacto moderado
    ka=1.25; 

    % Factor de tamaño según el modulo
    ks=2E-05*m^4 - 0.0011*m^3 + 0.0164*m^2 - 0.0647*m + 1.06;
    
    % Factor de distribucion de carga (montaje menos rígido)
    %km=5E-11*b^4 - 4E-08*b^3 + 1E-05*b^2 - 0.0007*b + 1.3074;
    km=3E-11*b^4 - 4E-08*b^3 + 1E-05*b^2 - 0.0007*b + 1.6071;

    kl=1;  %Factor engranaje loco
    kb=1;  %Factor espesor aro para engranaje solido
    
    %Factor dinamico    
    Qv=11;              %Indice de calidad del engranaje
    B=(12-Qv)^(2/3)/4;
    A=50+56*(1-B);
    
    %Vel tangencial [m/s]
    Vt1=(Rpm_max*2*pi/60)*R_p11/1000; %Par engranajes 1
    Vt2=(Rpm_max*2*pi/60)*R_p21/1000; %Par engranajes 2

    kv1=(A/(A+sqrt(200*Vt1)))^B;
    kv2=(A/(A+sqrt(200*Vt2)))^B;

    % Factor de Lewis
    Y1 = -0.00000005*z11^4+0.00001*z11^3-0.0007*z11^2+0.0241*z11+0.05;
    Y2 = -0.00000005*zmin^4+0.00001*zmin^3-0.0007*zmin^2+0.0241*zmin+0.05;
    
    %% Esfuerzo flexión con factores de comportamiento
    sigma_f1=abs(U1/(b*m*Y1)*(ka*ks*km*kl*kb/kv1)); %Esfuerzo en [N/mm^2 o MPa] 
    sigma_f2=abs(U2/(b*m*Y2)*(ka*ks*km*kl*kb/kv2)); %Esfuerzo en [N/mm^2 o Mpa] 

    % Asumiendo temperatura en caja de 80° C
    kt=(460+(80*9/5)+32)/620;   %Factor temperatura
    kr=1;                    %Factor fiabilidad    
    kn=9.4518*(10^6)^(-0.148);  %Factor duracion
    sigma_adm=St*kn/(kt*kr); %Esfuerzo admisible
   
        if (sigma_adm/sigma_f1<NS) || (sigma_adm/sigma_f2<NS)
            m=m+0.5;   
        end
    end
end

% Especificaciones de los engranajes
function Datos=Esp_engranajes(z, m)
    
%Parametros definidos, revisar función Calculo_engranajes
beta=23;    %Angulo de helice
alpha=20;   %Angulo de presión
psi=12;     %Factor guiado asumiendo calidad y condiciones normales

%Parametros calculados en [mm]

b=psi*m;  %Ancho de diente
p=pi*m;   %Paso del diente
j=0.25*m; %Holgura

h_c=m;     %Addendum
h_f=m+j;   %Deddendum
h=h_c+h_f; %Altura del diente

S=p/2; %Espesor

R_prim=(m/2)*(z/cosd(beta)); %Radio primitivo
R_c=R_prim+h_c;               %Radio cabeza
R_f=R_prim-h_f;              %Radio fondo
R_b=R_prim*cosd(alpha);      %Radio base

%d=m/2*zt/cosd(beta);       %Distancia entre ejes

Datos=[z, m, beta, alpha, b, p, R_prim, h_c, h_f, j, h, S, R_c, R_f, R_b];

end

% Diametro eje principal
function [d1e, R1_Ay, R1_Az, R2_Ay, R2_Az, R1_Bx, R1_By, R1_Bz, R2_Bx, R2_By, R2_Bz]=Eje_principal(U1, Fr1, Fa1, z11, U2, Fr2, Fa2, z21, m, S_ut)
    
    psi=12;
    beta=23;
    
    %Asumiendo distancias [mm]
    a=70;    %Distancia entre engranajes (2*b/2 ya incluida)
    b=m*psi; %Ancho engranaje
    L_r=20;  %Mitad ancho rodamiento+ buje (B/2 + buje)
    
    % Radios primitivos
    Rp1=(m/2)*(z11/cosd(beta)); %Radio primitivo primer marcha
    Rp2=(m/2)*(z21/cosd(beta)); %Radio primitivo segunda marcha
    
    %% Primera relacion, solo se considera ese par de engranajes
    %Balance fuerzas
    
    %Fuerzas en X
    R1_Bx=Fa1;
    
    %Momentos en Z respecto A
    R1_By=(Fr1*(b/2+L_r)+ Fa1*Rp1) / (2*(b/2+L_r) +a);
    
    %Momentos en Y respecto A
    R1_Bz=(U1*(b/2+L_r)) / (2*(b/2+L_r) +a);
     
    R1_Ay=Fr1-R1_By; %Fuerzas en Y
    R1_Az=U1-R1_Bz;  %Fuerzas en Z
    
    %Cortantes XY
    V11_xy=R1_Ay;        %Cortante A al engranaje
    V12_xy=V11_xy-Fr1;   %Cortante del engranaje a B
    V13_xy=V12_xy+R1_By; %Comprobacion (0)
    
    %Momentos Z
    M11_z=V11_xy*(b/2+L_r); 
    M12_z=M11_z+(Fa1*Rp1);            %Momento maximo
    M13_z=M12_z+V12_xy*((b/2+L_r)+a); %Comprobacion (0 prox)

    %Cortantes XZ
    V11_xz=R1_Az;        %Cortante A al engranaje
    V12_xz=V11_xz-U1;    %Cortante del engranaje a B
    V13_xz=V12_xz+R1_Bz; %Comprobacion (0 prox)
    
    %Momentos Y
    M11_y=V11_xz*(b/2+L_r);           %Momento maximo
    M12_y=M11_y+V12_xz*((b/2+L_r)+a); %Comprobacion (0 prox)
    
    %Torque y momento máximo
    M1=sqrt(M12_z^2 + M11_y^2); %Momento en el punto maximo
    T1=U1*Rp1;                  %Torque en el punto maximo

    %Momento en rodamientos
    %Rodamiento en A
    %M1_A=sqrt( (V11_xy*L_r)^2 + (V11_xz*L_r)^2);    
    %M1_A=0;
    %T1_A=-T1; %Torque en A
    %Rodamiento en B
    %M1_B=sqrt( (M11_y+V12_xz*(b/2+a))^2 + (M12_z+V12_xy*(b/2+a))^2);    
    %M1_B=0;
    %T1_B=0; %Torque en B

    %% Segunda relacion, solo se considera ese par de engranajes
    %Balance fuerzas

    %Fuerzas en X
    R2_Bx=Fa2;

    %Momentos en Z respecto A
    R2_By=(Fr2*(b/2+L_r+a)+ Fa2*Rp2) / (2*(b/2+L_r) +a);

    %Momentos en Y respecto A
    R2_Bz=(U2*(b/2+L_r+a)) / (2*(b/2+L_r) +a);
    
    R2_Ay=Fr2-R2_By;  %Fuerzas en Y
    R2_Az=U2-R2_Bz;   %Fuerzas en Z
    
    %Cortantes XY
    V21_xy=R2_Ay;        %Cortante A al engranaje
    V22_xy=V21_xy-Fr2;   %Cortante del engranaje a B
    V23_xy=V22_xy+R2_By; %Comprobacion (0)
    
    %Momentos Z
    M21_z=V21_xy*(b/2+L_r+a); 
    M22_z=M21_z+(Fa2*Rp2);        %Momento maximo
    M23_z=M22_z+V22_xy*(b/2+L_r); %Comprobacion (0 prox)
    
    %Cortantes XZ
    V21_xz=R2_Az;        %Cortante A al engranaje
    V22_xz=V21_xz-U2;    %Cortante del engranaje a B
    V23_xz=V22_xz+R2_Bz; %Comprobacion (0)
    
    %Momentos Y
    M21_y=V21_xz*(b/2+L_r+a);     %Momento maximo
    M22_y=M21_y+V22_xz*(b/2+L_r); %Comprobacion (0 prox)
    
    %Torque y momentos en el punto máximo
    M2=sqrt(M22_z^2 + M21_y^2); %Momento en el punto maximo
    T2=-U2*Rp2;                 %Torque en el punto maximo

    %Momento en rodamientos
    %Rodamiento en A
    %M2_A=sqrt( (V21_xy*L_r)^2 + (V21_xz*L_r)^2);   
    %M2_A=0;
    %T2_A=T2; %Torque en A
    %Rodamiento en B
    %M2_B=sqrt( (M21_y+V22_xz*b/2)^2 + (M22_z+V22_xy*b/2)^2);    
    %M2_B=0;
    %T2_B=0; %Torque en B

    %% Diametros

    NS=2; %Factor de seguridad

    %Factores para esfuerzo a flexión
    k_t=2;    %Asumiendo concentrador de esfuerzos
    q_f=0.78; %Asumiendo un radio de muesca de 1.5mm
    k_fm=1+q_f*(k_t-1);
    
    %Factores para esfuerzo a torsión
    q_t=0.8;
    k_fsm=1+q_t*(k_t-1);

    %Factor de corrección para el esfuerzo a la fatiga
    d=20; %Diametro asumido inicial
    d_1m=0;  d_2m=0; 

    C_carga=1;    %Carga a flexión
    C_sup=0.77;   %Mecanizado
    C_temp=1;     %Temperatura menor a 450 °C
    C_conf=0.868; %95% de confiabilidad

    Se=S_ut/2;  
    
    %Se asume un diametro, se obtiene el factor de tamaño, se obtiene el diametro segun esfuerzos y se corrige el diametro asumido
    %El ciclo se repite hasta que el diametro del factor es igual al obtenido por esfuerzos
    while roundn(d,-1)~=roundn(d_1m,-1) 
    C_tam=1.189*d^(-0.097);
    S_e=C_carga*C_tam*C_sup*C_temp*C_conf*Se;
    
    d_1m=( 32*NS/(pi) * ( (k_fm*M1)/S_e + sqrt((3/4)*(k_fsm*T1)^2)/S_ut ) )^(1/3);
    if ( roundn(d_1m,-1)~= roundn(d,-1))
        d=d_1m; d_1m=0;
    end
    end
    
    while roundn(d,-1)~=roundn(d_2m,-1) 
    C_tam=1.189*d^(-0.097);
    S_e=C_carga*C_tam*C_sup*C_temp*C_conf*Se;

    d_2m=( 32*NS/(pi) * ( (k_fm*M2)/S_e + sqrt((3/4)*(k_fsm*T2)^2)/S_ut ) )^(1/3);
    if (roundn(d_2m,-1)~= roundn(d,-1))
        d=d_2m;    d_2m=0;        
    end
    end

    %Comparando diametro mayor y normalizando el eje a entero
    if (d_1m>d_2m) % Se escoge el diametro mayor - se aproxima entero mayor
        d1e=ceil(d_1m);
    else
        d1e=ceil(d_2m);
    end
    
end

% Diametro eje secundario
function [d2e, R1_Ax, R1_Ay, R1_Az, R2_Ax, R2_Ay, R2_Az, R1_By, R1_Bz, R2_By, R2_Bz]=Eje_secundario(U1, Fr1, Fa1, z12, U2, Fr2, Fa2, z22, m, S_ut)
 
    %Asumiendo distancias [mm]
    a=70;   %Distancia entre engranajes (2*b/2 ya incluida)
    b=25;   %Ancho engranaje
    L_r=20; %Mitad ancho rodamiento+ buje (B/2 + buje)

    % Radios primitivos
    beta=23;
    Rp1=(m/2)*(z12/cosd(beta)); %Radio primitivo primer marcha
    Rp2=(m/2)*(z22/cosd(beta)); %Radio primitivo segunda marcha
  
    %% Primera relacion, solo se considera ese par de engranajes
    %Balance fuerzas

    %Fuerzas en X
    R1_Ax=Fa1;

    %Momentos en Z respecto A
    R1_By=(Fr1*(b/2+L_r)-Fa1*Rp1) / (2*(b/2+L_r) +a);

    %Momentos en Y respecto A
    R1_Bz=(U1*(b/2+L_r)) / (2*(b/2+L_r) +a);
    
    R1_Ay=Fr1-R1_By;  %Fuerzas en Y
    R1_Az=U1-R1_Bz;   %Fuerzas en Z
    
    %Cortantes XY
    V11_xy=-R1_Ay;       %Cortante A al engranaje
    V12_xy=V11_xy+Fr1;   %Cortante del engranaje a B
    V13_xy=V12_xy-R1_By; %Comprobacion (0)
    
    %Momentos Z
    M11_z=V11_xy*(b/2+L_r);           %Momento maximo 
    M12_z=M11_z+(Fa1*Rp1); 
    M13_z=M12_z+V12_xy*((b/2+L_r)+a); %Comprobacion (0 prox)

    %Cortantes XZ
    V11_xz=-R1_Az;       %Cortante A al engranaje
    V12_xz=V11_xz+U1;    %Cortante del engranaje a B
    V13_xz=V12_xz-R1_Bz; %Comprobacion (0 prox)
    
    %Momentos Y
    M11_y=V11_xz*(b/2+L_r);           %Momento maximo
    M12_y=M11_y+V12_xz*((b/2+L_r)+a); %Comprobacion (0 prox)
    
    %Torque y momentos en el punto máximo
    M1=sqrt(M11_z^2 + M12_y^2); %Momento en el punto máximo
    T1=U1*Rp1;                  %Torque en el punto máximo

    %Momento en rodamientos
    %Rodamiento en A
    %M1_A=sqrt( (V11_xy*L_r)^2 + (V11_xz*L_r)^2); 
    %M1_A=0;
    %T1_A=0; %Torque en A
    %Rodamiento en B
    %M1_B=sqrt( (M11_y+V12_xz*(b/2+a))^2 + (M12_z+V12_xy*(b/2+a))^2);   
    %M1_B=0;
    %T1_B=T1; %Torque en B

    %% Segunda relacion, solo se considera ese par de engranajes
    %Balance fuerzas

    %Fuerzas en X
    R2_Ax=Fa2;

    %Momentos en Z respecto A
    R2_By=(Fr2*(b/2+L_r+a)-Fa2*Rp2) / (2*(b/2+L_r) +a);

    %Momentos en Y respecto A
    R2_Bz=(U2*(b/2+L_r+a)) / (2*(b/2+L_r) +a);
    
    R2_Ay=Fr2-R2_By; %Fuerzas en Y
    R2_Az=U2-R2_Bz;  %Fuerzas en Z
    
    %Cortantes XY
    V21_xy=-R2_Ay;       %Cortante A al engranaje
    V22_xy=V21_xy+Fr2;   %Cortante del engranaje a B
    V23_xy=V22_xy-R2_By; %Comprobacion (0)
    
    %Momentos Z
    M21_z=V21_xy*(b/2+L_r+a);     %Momento maximo
    M22_z=M21_z+(Fa2*Rp2); 
    M23_z=M22_z+V22_xy*(b/2+L_r); %Comprobacion (0 prox)
    
    %Cortantes XZ
    V21_xz=-R2_Az;       %Cortante A al engranaje
    V22_xz=V21_xz+U2;    %Cortante del engranaje a B
    V23_xz=V22_xz-R2_Bz; %Comprobacion (0)
    
    %Momentos Y
    M21_y=V21_xz*(b/2+L_r+a);     %Momento maximo
    M22_y=M21_y+V22_xz*(b/2+L_r); %Comprobacion (0 prox)
    
    %Torque y momentos en el punto máximo
    M2=sqrt(M22_z^2 + M21_y^2); %Momento en el punto máximo
    T2=U2*Rp2;                  %Torque en el punto máximo

    %Momento en rodamientos
    %Rodamiento en A
    %M2_A=sqrt( (V21_xy*L_r)^2 + (V21_xz*L_r)^2);   
    %M2_A=0;
    %T2_A=0; %Torque en A
    %Rodamiento en B
    %M2_B=sqrt( (M21_y+V22_xz*b/2)^2 + (M22_z+V22_xy*b/2)^2); 
    %M2_B=0;
    %T2_B=T2; %Torque en B

    %% Diametros

    NS=2; %Factor de seguridad

    %Factores para esfuerzo a flexión
    k_t=2;    %Asumiendo concentrador de esfuerzos
    q_f=0.78; %Asumiendo un radio de muesca de 1.5mm
    k_fm=1+q_f*(k_t-1);
    
    %Factores para esfuerzo a torsión
    q_t=0.8;
    k_fsm=1+q_t*(k_t-1);

    %Factor de corrección para el esfuerzo a la fatiga

    d=20; %Diametro asumido inicial
    d_1m=0;  d_2m=0; 

    C_carga=1;    %Carga a flexión
    C_sup=0.77;   %Mecanizado
    C_temp=1;     %Temperatura menor a 450 °C
    C_conf=0.868; %95% de confiabilidad

    Se=S_ut/2;  
    
    %Se asume un diametro, se obtiene el factor de tamaño, se obtiene el diametro segun esfuerzos y se corrige el diametro asumido
    %El ciclo se repite hasta que el diametro del factor es igual al obtenido por esfuerzos
    while roundn(d,-1)~=roundn(d_1m,-1) 
    C_tam=1.189*d^(-0.097);
    S_e=C_carga*C_tam*C_sup*C_temp*C_conf*Se;
    
    d_1m=( 32*NS/(pi) * ( (k_fm*M1)/S_e + sqrt((3/4)*(k_fsm*T1)^2)/S_ut ) )^(1/3);
    if ( roundn(d_1m,-1)~= roundn(d,-1))
        d=d_1m; d_1m=0;
    end
    end
    
    while roundn(d,-1)~=roundn(d_2m,-1) 
    C_tam=1.189*d^(-0.097);
    S_e=C_carga*C_tam*C_sup*C_temp*C_conf*Se;

    d_2m=( 32*NS/(pi) * ( (k_fm*M2)/S_e + sqrt((3/4)*(k_fsm*T2)^2)/S_ut ) )^(1/3);
    if (roundn(d_2m,-1)~= roundn(d,-1))
        d=d_2m; d_2m=0;
    end
    end
    
    %Comparando diametro mayor y normalizando el eje a entero
    if (d_1m>d_2m) %Se escoge el diametro mayor - se aproxima entero mayor
        d2e=ceil(d_1m);
    else
        d2e=ceil(d_2m);
    end

end

% Embrague de disco
% Se debe normalizar el resultado (verificar estriado)
function embrague=Calculo_embrague(Motores, Numero_motor)     
        
            N=1.5; %Factor de seguridad vehiculo cargas normales
            Par_d=Motores(Numero_motor,4)*N*(100/9.81);
        
            P_max=2.5; %Se asume presión vehiculos turismo
            mu=0.4;    %Se asume coeficiente de rozamiento
            R_ext=(Par_d/(2.75*P_max*mu))^(1/3);
            R_int=0.7*R_ext;
        
            %Calculo embrague desgastado
            %P=P_max*R_int/R_ext;
            F=2*pi*P_max*R_int*(R_ext-R_int);
        
            c=2; %Caras en contacto
            T_roz=c*mu*F/2*(R_ext+R_int)*(9.81/100);
        
            if(T_roz>Motores(Numero_motor,4)) 
                D_min=2*R_ext*10; %Diametro minimo embrague en [mm]
            else
                D_min=0; %El embrague no es capaz de transmitir el par completo
            end
            
            embrague=D_min;
end

% Rodamientos
function rodamiento_b=calc_rod_bolas(R1_y, R1_z, R2_y, R2_z, rpm)
    
    a=3; %Rodamientos bolas

    T_vida=250000; %Tiempo de vida caja estimado [km]
    Vel_media=60;  %Velocidad media del vehiculo [km/h]
    T_duracion=T_vida/Vel_media; %Horas duración caja

    F_1r= sqrt(R1_y^2 + R1_z^2);
    F_2r= sqrt(R2_y^2 + R2_z^2);

    if(F_1r>=F_2r)
        P_radial=F_1r;
    else
        P_radial=F_2r;
    end

    L_10=((T_duracion*rpm*60)/10^6) / (0.02+4.439*log(1/0.9)^(1/1.483));
    C_radial= (L_10)^(1/a) *1.2*P_radial;

    rodamiento_b=C_radial;
        
end

function rodamiento_r=calc_rod_rodillos(R1_x, R1_y, R1_z, R2_x, R2_y, R2_z, rpm)

load('VidaUtil.mat') %Se carga el archivo de vida util engranajes segun porcentaje

a=10/3; %Rodamientos de rodillos
%e=0.2; %para rodamientos de las series 10,18, 19, 2, 3 y 4
e=0.3; %para rodamientos de las series 12, 20, 22, 23, 28, 29, 30 y 39

    T_vida=250000; %Tiempo de vida caja estimado [km]
    Vel_media=60;  %Velocidad media del vehiculo [km/h]
    T_duracion=T_vida/Vel_media; %Horas duración caja    

    F_1cr= sqrt(R1_y^2 + R1_z^2);
    F_2cr= sqrt(R2_y^2 + R2_z^2);

    if (abs(R1_x)/F_1cr) <e
        P_1=(0.4*R1_x+0.92*F_1cr);
    else
        P_1=F_1cr;
    end

    if (abs(R2_x)/F_2cr) <e
        P_2=(0.4*R2_x+0.92*F_2cr);
    else
        P_2=F_2cr;
    end

    P_combinada=(P_1^3*Vida_util(1,2) + P_2^3*Vida_util(2,2))^(1/3);
    L_10=((T_duracion*rpm*60)/10^6) / (0.02+4.439*log(1/0.9)^(1/1.483));
    C_combinada= (L_10)^(1/a) *1.2*P_combinada;

    rodamiento_r=C_combinada;
end

