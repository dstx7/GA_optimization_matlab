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
load('InfoMotoresDastan.mat');

%%==============================================%%

[R_rodadura, R_pendiente, R_aerodinamica, R_inercia, P_min]=Calculo_Pot_min(V, W, a);

%%Estructura resultado
% || ID || Primera || Segunda || Rel Reductor || Rel Diferencial
Resultados=Calculo_Relaciones(Motores, V, rueda_carga, W, R_rodadura, R_pendiente, R_aerodinamica, R_inercia);
% || ID || Primera || Segunda || Rel Reductor || Rel Diferencial || Dmin embrague
Resultados=[Resultados, Calculo_embrague(Motores)]
writematrix(Resultados, 'Resultados.xlsx');

%% Funciones

%Resistencias y potencia mínima [KW]
function  [Fr, Fp, Fa, Fj, Pot_min]=Calculo_Pot_min(Vel_max, Peso, Aceleracion)
    
    V = Vel_max/3.6;
    
    %%resistencia a la rodadura
    fo = 0.009;  
    fs = 0.004;
    Cfr = fo +  3.24*fs*(V/(2.237*100))^2.5; %coeficiente de resistencia a la rodadura 
    Fr = Cfr*Peso;     
    
    %%resistencia aerodinamica
    CD = 0.39;      %coeficiente de resistencia aerodinamica
    A = 1.40875;    %area frontal del vehiculo 
    rho_a = 1.204;  %densidad del aire a T ambiente
    Fa = CD*A*(rho_a/2)*V^2;  

    %%resistencia pendiente
    pendiente=0.12; %aproximación de pendiente máxima
    Fp=Peso*pendiente;

    %%resistencia inercia
    j=Aceleracion;
    Fj=(Peso/9.81)*j;
    
    %%Potencia mínima
    Pot_min = (Fa + Fr)*V/1000; %Potencia resultante sin pendiente
end

%Relaciones de la caja
function cambios=Calculo_Relaciones(Motores, Vel_max, R_rueda_carga, Peso, Fr, Fp, Fa, Fj)

Dimension_matriz=size(Motores);  %Motores es el archivo "InfoMotores"
cambios=[0, 0, 0, 0, 0];

%Constantes 
N_eficiencia=0.96; %Factor seguridad por perdidas de eficiencia
miu_coeff = 0.8; %coeficiente de friccion estatica en condiciones de superficie seca a 90m/h
alpha_st = 16.69; %Angulo máximo de inclinación de la carretera
V_subiendo = 15; %velocidad en ascenso a pendiente de 30
rho_air = 1.204; % densidad del aire a 20 °C [kg/m^3]
CD = 0.39; % coeficiente de resistencia aerodinamica
A= 1.40875; %area frontal del vehiculo [m^2]
rho_air = 1.204; % densidad del aire a 20 °C [kg/m^3]

n_r=(Vel_max*60/3.6)/(pi*R_rueda_carga*2);

%Test primera marcha
R_resistente_primera=Fr+Fp+Fj; 
T_resistente_primera=R_resistente_primera*R_rueda_carga; 

%Test segunda marcha
R_resistente_segunda=Fr+Fa;
T_resistente_segunda=R_resistente_segunda*R_rueda_carga;

I_E=5.13; %Constante relación transmisión final optimizar torque

    for i= 1:Dimension_matriz(1)
        relacion_reductor=Vel_max/(Motores(i,2)*2*pi*R_rueda_carga*60*10^-3);

        I_primera_max=(R_rueda_carga*Peso*miu_coeff)/(Motores(i,4)*N_eficiencia); %Se debe dividir entre I_E?
        %I_primera_min= (R_rueda_carga*Peso*((Fr/Peso)*cosd(alpha_st)+sind(alpha_st))+0.5*CD*rho_air*A*V_subiendo^2)/(Motores(i,4)*N_eficiencia);
        I_primera_min= (R_rueda_carga*Peso*((Fr/Peso)*cosd(alpha_st)+sind(alpha_st)))/(Motores(i,4)*N_eficiencia*I_E);
        
        I_segunda_max= (3.6*pi*R_rueda_carga*Motores(i,2))/(30*Vel_max*I_E);
        I_segunda_min= (3.6*pi*R_rueda_carga*Motores(i,3))/(30*Vel_max*I_E);
        
        marcha=Opt_cambios(Motores(i,2),Motores(i,4), N_eficiencia, n_r, I_primera_min,I_primera_max, ...
            I_segunda_min, I_segunda_max , T_resistente_primera, T_resistente_segunda);
        
        relacion_diferencial=Motores(i,2)/(n_r*marcha.y);
        cambios=[cambios; Motores(i,1), marcha.x, marcha.y, relacion_reductor, relacion_diferencial];
    end
    
end

%Embrague de disco
%Se debe normalizar el resultado (verificar estriado)
function embrague=Calculo_embrague(Motores)     

embrague=[0];
Dimension_matriz=size(Motores);

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

        embrague=[embrague; D_min];
    end
end

%Optimización cambios
function Obj_cambio=Opt_cambios(Rpm_max, T_max, Factor_seguridad, N_r, Rel_primera_min,Rel_primera_max, Rel_segunda_min, Rel_segunda_max , Primera_test, Segunda_test)
   

% Create optimization variables

x4 = optimvar("x","LowerBound",Rel_primera_min); %Solo limite inferior
%Si se desea con limite superior y inferior
%x4 = optimvar("x","LowerBound",Rel_primera_min,"UpperBound",Rel_primera_max); 
y4 = optimvar("y","LowerBound",Rel_segunda_min,"UpperBound",Rel_segunda_max);

% Set initial starting point for the solver
initialPoint3.x = (Rel_primera_min+Rel_primera_max)/2;
initialPoint3.y = (Rel_segunda_max+Rel_segunda_min)/2;

% Create problem
problem = optimproblem("ObjectiveSense","Maximize");

% Define problem objective

rel_diferencial=Rpm_max/(N_r*y4);
T_rueda_primera=(T_max*Factor_seguridad)/(x4/rel_diferencial);
T_rueda_segunda=(T_max*Factor_seguridad)/(y4/rel_diferencial);
problem.Objective = (T_rueda_primera-Primera_test)+(T_rueda_segunda-Segunda_test);

% Display problem information
show(problem);

% Solve problem
[solution,objectiveValue,reasonSolverStopped] = solve(problem,initialPoint3);

% Clear variables
clearvars x4 y4 initialPoint3 reasonSolverStopped objectiveValue;

% Display results
Obj_cambio=solution;

end


