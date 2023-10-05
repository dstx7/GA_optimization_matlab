%% Datos fijos

%Llanta 125/80 R13
f=0.9; %Factor de compensación llanta bajo carga
%rueda_carga=f*(13*25.4 +0.8*125*2)/(2*1000) %Radio de la llanta bajo carga
rueda_carga=0.265;

a=50*(1/(3.6*8)); %aceleración: 0 a 50 [km/h] en 8 segundos
V=80; % Velocidad máxima [km/h]
W = 700*9.81;  %Peso total del vehiculo en kg (antes tenia 750 se modifico por los objetivos)

% Cargar tabla información motores
% Estructura del archivo Infomotores
% || Id motor || Rpm max || Rpm min || Par max ||
load('InfoMotores.mat');

%% ============================================== %%

% Calculo de resistencias y potencia mínima
[R_rodadura, R_pendiente, R_aerodinamica, R_inercia, P_min]=Calculo_Pot_min(V, W, a);

% Estructura resultado
% || ID || Primera min || Primera_max || Segunda min || Segunda max
Intervalos=Calculo_Relaciones(Motores, V, rueda_carga, W, R_rodadura);

% Resultados piezas
% || ID || Dmin embrague
Componentes=[Calculo_embrague(Motores)];
% writematrix(Resultados, 'Resultados.xlsx');


%% Funciones

% Resistencias y potencia mínima [KW]
function  [Fr, Fp, Fa, Fj, Pot_min]=Calculo_Pot_min(Vel_max, Peso, Aceleracion)
    
    V = Vel_max/3.6;
    
    %% resistencia a la rodadura
    fo = 0.009;  
    fs = 0.004;
    Cfr = fo +  3.24*fs*(V/(2.237*100))^2.5; %coeficiente de resistencia a la rodadura 
    Fr = Cfr*Peso;     
    
    %% resistencia aerodinamica
    CD = 0.39;      %coeficiente de resistencia aerodinamica
    A = 1.40875;    %area frontal del vehiculo 
    rho_a = 1.204;  %densidad del aire a T ambiente
    Fa = CD*A*(rho_a/2)*V^2;  

    %% resistencia pendiente
    pendiente=0.12; %aproximación de pendiente máxima
    Fp=Peso*pendiente;

    %% resistencia inercia
    j=Aceleracion;
    Fj=(Peso/9.81)*j;
    
    %% Potencia mínima
    Pot_min = (Fa + Fr)*V/1000; %Potencia resultante sin pendiente
end

% Intervalos relaciones de la caja
function cambios=Calculo_Relaciones(Motores, Vel_max, R_rueda_carga, Peso, Fr)

    Dimension_matriz=size(Motores);  %Motores es el archivo "InfoMotores"
    cambios=zeros(Dimension_matriz(1),5);
    
    %% Constantes
    N_eficiencia=0.96; %Factor seguridad por perdidas de eficiencia
    miu_coeff = 0.8; %coeficiente de friccion estatica en condiciones de superficie seca a 90m/h
    alpha_st = 16.69; %Angulo máximo de inclinación de la carretera
    
    I_E=3.73; %Constante relación transmisión final-optimizar torque

    for i= 1:Dimension_matriz(1)
        
        %Ecuaciones de: 
        % "Design and Optimization of a Drivetrain with Two-speed Transmission for Electric Delivery Step Van"
        I_primera_min= (R_rueda_carga*Peso*((Fr/Peso)*cosd(alpha_st)+sind(alpha_st)))/(Motores(i,4)*N_eficiencia*I_E);
                        
        I_segunda_max= (3.6*pi*R_rueda_carga*Motores(i,2))/(30*Vel_max*I_E);
        I_segunda_min= (3.6*pi*R_rueda_carga*Motores(i,3))/(30*Vel_max*I_E);
        
        %Ecuacion de:
        % "Comprehensive design and optimization of an electric vehicle powertrain equipped with a two-speed dual-clutch transmission"
        I_primera_max=(R_rueda_carga*Peso*miu_coeff)/(Motores(i,4)*N_eficiencia); %Se debe dividir por I_E?
        
        cambios(i,1)=Motores(i,1); cambios(i,2)=I_primera_min; cambios(i,3)=I_primera_max;
        cambios(i,4)=I_segunda_min; cambios(i,5)=I_segunda_max;
    end
    
end

% Engranajes (dientes,módulo y distancia entre ejes)
function engranajes=Calculo_engranajes(I_primera, I_segunda, Potencia, Rpm_max, Tor_max)
    
    %% Número de dientes y relaciones reales
    % Consideraciones:
    % Todos engranajes misma distancia entre centros
    % Engranajes con el mismo módulo
    % Constante la suma de número de dientes de cada par de engranajes
    % Las relaciones se ingresan en valores >1 
    
    beta=20;    %Angulo de helice
    alpha=20;   %Angulo de presión
    zn=14;       %Zn >=14 número de dientes virtual
    
    %% Constantes del material   " ¿20MnCr5? "
    k_adm=80; %Acero aleado, cementado y templado [kg/cm^2]
    sigma_adm=3400; %Esfuerzo admisible [kg/cm^2]
    %% 

    %Primer conjunto
    z11=ceil(zn*cosd(beta)^3); %Número de dientes real - se escoge el mínimo
    z12=round(z11*I_primera);
    cons=z11+z12;
    i1=z12/z11;  %Primera relación dada por los números de dientes
    
    %Segundo conjunto
    z21=round(cons/(1+(1/I_segunda)));
    z22=cons-z21;
    i2=z21/z22; %Segunda relación dada por los números de dientes
    
    %% Módulo y distancia entre ejes
    
    T_vida=250000; % Tiempo de vida caja estimado [km]
    Vel_media=60; %Velocidad media del vehiculo [km/h]
    
    T_duracion=T_vida/Vel_media; %Horas duración caja
    
    load('VidaUtil.mat') %Se carga el archivo de vida util engranajes segun porcentaje
    Dimension_matriz=size(Vida_util);
        
    for i= 1:Dimension_matriz(1)
        Vida_util(i,3)=T_duracion*Vida_util(i,2); %Horas utiles de cada marcha
        phi=16.895*Vida_util(i,3)^(-0.332);
        Vida_util(i,4)=phi*k_adm;
    end

    psi=10; %Factor guiado asumiendo calidad y condiciones normales
    Pot=Potencia*1.34102

    % Primer conjunto - modulo en [mm]
    m1=(143240*Pot*(i1+1)*(cosd(20)^4)/(Vida_util(1,4)*z11^2*psi*cosd(20)*sind(20)*Rpm_max*i1))^(1/3)*10;

    % Segundo conjunto - modulo en [mm]
    % El modulo se halla a partir de los dientes del piñon
    if (z21<z22)
        zmin=z21; %El modulo se halla a partir de los dientes del piñon
    else
        zmin=z22;
    end
    m2=(143240*Pot*(i2+1)*(cosd(20)^4)/(Vida_util(2,4)*zmin^2*psi*cosd(20)*sind(20)*Rpm_max*i2))^(1/3)*10;
    
    if (m1>m2) %Se escoge el modulo mayor - se aproxima entero mayor
        m=ceil(m1);
    else
        m=ceil(m2);
    end
    
    b=m*psi;
    d=m/2*cons/cosd(20); %distancia entre ejes

    %% Fuerzas y comprobación sobre los engranajes 
    % Esfuerzo a flexión
    
    T_max=Tor_max*1000; %Torque maximo [Nmm]

    % Primera marcha
    Rp1=m/2*(z11/cosd(beta)); %Radio primitivo piñon
    U1=T_max/Rp1;
    Fr1=U1*tand(alpha);
    Fa1=U1*tand(beta);
    W1=sqrt(U1^2+Fr1^2+Fa1^2);
    
    Y1 = -5E-11*z11^6 + 1E-08*z11^5 - 2E-06*z11^4 + 0.0001*z11^3 - 0.0039*z11^2 + 0.0706*z11 - 0.2013;
    sigma_f1=U1/(b*m*Y1)*10.197; %Esfuerzo de [N/mm^2] a [kgf/cm^2]

    % Segunda marcha
    Rp2=m/2*(zmin/cosd(beta)); %Radio primitivo piñon  
    U2=T_max/Rp2;
    Fr2=U2*tand(alpha);
    Fa2=U2*tand(beta);
    W2=sqrt(U2^2+Fr2^2+Fa2^2);
    
    Y2 = -5E-11*zmin^6 + 1E-08*zmin^5 - 2E-06*zmin^4 + 0.0001*zmin^3 - 0.0039*zmin^2 + 0.0706*zmin - 0.2013;
    sigma_f2=U2/(b*m*Y2)*10.197; %Esfuerzo de [N/mm^2] a [kgf/cm^2]
    
    %Se evalua si los esfuerzos en los engranajes son mayores a los admisibles
    if (sigma_f1>sigma_adm || sigma_f2>sigma_adm) 

        m=999999; d=1; beta=1;
        z11=1; z12=1; z21=1; z22=1;

    end

    engranajes=[ z11, z12, z21, z22, m, d, beta];
end

% Engranajes (dimensiones)
function dim_engranajes=Calculo_dimensiones_engranajes( z11, z12, z21, z22, m, d, beta)
    
    psi=10;
    
    b=psi*m;

end

% Embrague de disco
% Se debe normalizar el resultado (verificar estriado)
function embrague=Calculo_embrague(Motores)     

    Dimension_matriz=size(Motores);
    embrague=zeros(Dimension_matriz(1),2);
    
        for i= 1:Dimension_matriz(1)
    
            N=1.5; %Factor de seguridad vehiculo cargas normales
            Par_d=Motores(i,4)*N*(100/9.81);
        
            P_max=2.5; %Se asume presión vehiculos turismo
            mu=0.4; %Se asume coeficiente de rozamiento
            R_ext=(Par_d/(2.75*P_max*mu))^(1/3);
            R_int=0.7*R_ext;
        
            %Calculo embrague desgastado
            P=P_max*R_int/R_ext;
            F=2*pi*P_max*R_int*(R_ext-R_int);
        
            c=2; %Caras en contacto
            T_roz=c*mu*F/2*(R_ext+R_int)*(9.81/100);
        
            if(T_roz>Motores(i,4)) 
                D_min=2*R_ext*10; %Diametro minimo embrague en [mm]
            else
                D_min=0; %El embrague no es capaz de transmitir el par completo
            end
            
            embrague(i,1)=Motores(i,1); embrague(i,2)=D_min;
        end
end

% Optimización cambios
function Obj_cambio=Opt_cambios(Rpm_max, T_max, Factor_seguridad, N_r, Rel_primera_min,Rel_primera_max, Rel_segunda_min, Rel_segunda_max)
    
    %% Con solo limite inferior en la primera marcha
    
    x4 = optimvar("x","LowerBound",Rel_primera_min); %Solo limite inferior
    y4 = optimvar("y","LowerBound",Rel_segunda_min,"UpperBound",Rel_segunda_max);
    
    % Set initial starting point for the solver
    initialPoint3.x = (Rel_primera_min+0.5)/2;
    initialPoint3.y = (Rel_segunda_max+Rel_segunda_min)/2;
    
    %% Con limite inferior y superior en la primera marcha
    
    %x4 = optimvar("x","LowerBound",Rel_primera_min,"UpperBound",Rel_primera_max);  % Si se desea con limite superior y inferior
    %y4 = optimvar("y","LowerBound",Rel_segunda_min,"UpperBound",Rel_segunda_max);
    
    % Set initial starting point for the solver
    %initialPoint3.x = (Rel_primera_min+Rel_primera_max)/2;
    %initialPoint3.y = (Rel_segunda_max+Rel_segunda_min)/2;
    
    %% Create problem
    problem = optimproblem("ObjectiveSense","Maximize");
    
    % Define problem objective
    
    rel_diferencial=Rpm_max/(N_r*y4);
    T_rueda_primera=(T_max*Factor_seguridad)/(x4/rel_diferencial);
    T_rueda_segunda=(T_max*Factor_seguridad)/(y4/rel_diferencial);
    problem.Objective = 0.6*(T_rueda_primera+T_rueda_segunda);
    
    % Display problem information
    show(problem);
    
    % Solve problem
    [solution,objectiveValue,reasonSolverStopped] = solve(problem,initialPoint3);
    
    % Clear variables
    clearvars x4 y4 initialPoint3 reasonSolverStopped objectiveValue;
    
    % Display results
    Obj_cambio=solution;

end


