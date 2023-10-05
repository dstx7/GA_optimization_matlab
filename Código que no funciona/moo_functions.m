%%Funciones objetivo
function f= moo_functions (x)
    
    %Se renombran las variables a iterar
    I_primera=x(1); I_segunda=x(2);
    
    global Tor_max z11 z12 z21 z22 m rho_eng rho_eje d1e d2e Rpm_max Fr1
    
    v= 150;   %viscosidad cinematica del aceite a temperatura de sumidero en centistokes (SAE 80W90 A 40 °C)
    alpha=20; %Angulo de presión
    beta=23;  %Angulo helice
    psi=12;   %Factor guiado asumiendo calidad y condiciones normales (revisar función Calculo_engranajes)
    b=psi*m;  %Ancho de diente
    g=9.81;   %Gravedad

    I_E=3.77; %Constante relación transmisión final-optimizar torque (Ver General)

    a=70;   %Distancia entre engranajes (2*b/2 ya incluida)
    L_r=20; %Mitad ancho rodamiento + buje (B/2 + buje)
    
    %% Funcion objetivo minimizar peso
    V_eje=( (2*L_r*(pi*d1e^2)/4) + a*(pi*(d1e+5)^2)/4 + (2*b*(pi*d1e^2)/4) )+ ( (2*L_r*(pi*d2e^2)/4)+ a*(pi*(d2e+5)^2)/4+ (2*b*(pi*d2e^2)/4));
    V_engranajes=(pi/4 * (m^2/(cosd(beta))^2) *b* (z11^2+z12^2+z21^2+z22^2));
    f(1)= (rho_eng*V_engranajes+rho_eje*V_eje)/(1E9)*g; %[N]
    
    %% Funcion objetivo torque
    %Para maximizar multiplicamos por menos toda la funcion
    Factor_seguridad=0.95;

    T_rueda_primera=(Tor_max*Factor_seguridad)/(1/(I_E*I_primera));
    T_rueda_segunda=(Tor_max*Factor_seguridad)/(1/(I_E*I_segunda));

    f(2)= -T_rueda_segunda - T_rueda_primera ;

    %% Funcion objetivo perdidas de potencia
    z31=12; %Asumiendo dientes del diferencial

    Pm1 = perdidasEngranadoPm(z11, z12, m, Rpm_max, v, Tor_max, alpha, beta);
    Pm2 = perdidasEngranadoPm(z21,z22, m, Rpm_max, v, Tor_max, alpha, beta);
    Pm3 = perdidasEngranadoPm(z31, I_E, m, Rpm_max/I_primera, v, Tor_max*I_primera,alpha, beta);
    T_GM_1st = 9550*I_primera*Pm1/Rpm_max; %Nm
    T_GM_2st = 9550*I_segunda*Pm2/Rpm_max; %Nm
    T_GM_3st = 9550*I_E*Pm3/(Rpm_max/I_primera); %Nm

    Tb1 = bearingsPb(Fr1, d1e, Rpm_max); %Perdidas por los rodamientos
    Tb2 = bearingsPb(Fr1, d2e, Rpm_max); %Perdidas por los rodamientos

    T_1st_a = (Tor_max -2*Tb1) * I_primera ;
    T_1st_b = (T_1st_a - T_GM_1st - T_GM_2st*I_segunda - 2*Tb2)*I_E;
    T_salida = T_1st_b - T_GM_3st - Tb2;

    %f(2)=-T_salida;

end

%% Funciones utilizadas para obtener las funciones objetivo

function [M1] = bearingsPb(fr, di, n)
    %Perdidas en los rodamientos asumiendo rodamientos conicos simples 
    %fr y fa corresponden a cargar radial y axial (estoy asumiendo solo carga
    %radial mientras encontramos diametros de los engranajes)

    %Para un rodamiento skf de 35mm
    do = di*1.8;
    
    f1 = 0.004;       %Coeficiente de friccion tabla 2 (norma ISO)
    p1 = fr/2;        %Carga dinamica del rodamiento Newton
    dm = (di + do)/2; %Diametro interior y exterior del rodamiento en mm
    
    M1 = f1*(p1)*(dm)/1000; %Torque dependiente de la carga en N-m
    
    Pb = (M1*n)/9549; %Perdidas de potencia de los rodamientos en kW

end

function [Pm] = perdidasEngranadoPm(z1,z2,modulo,n1,v,Tm,alpha,betha)
    
    psi=12; %Factor guiado asumiendo calidad y condiciones normales (constante definida tambien en la funcion calculo engranajes)

    u=z2/z1;                       %Relación de la marcha
    modulo_c = modulo/cosd(betha); %Modulo aparente
    pd = 1/modulo;                 %Paso diametral mm^-1
    bw=modulo*psi;                 %Ancho de diente

    dp_1 = z1*modulo_c; %Diametro de paso piñon mm
    dp_2 = z2*modulo_c; %Diameto de paso engrane mm
    
    addn = 1/pd; %Addendum
    
    ro_2 = dp_2/2 + addn; %Radio exterior del engranaje (conducido) en mm
    rw_2 = dp_2/2;        %Radio de paso operativo en mm
    ro_1 = dp_1/2 + addn; %Radio exterior del piñon (conductor) en mm 
    rw_1 = dp_1/2;        %Radio de paso operativo en mm
    
    V_pitch = n1*(2*pi/60)*(rw_1/1000);
    
    C1 = 3.239; %Constantes de la norma
    j = -0.223;
    g = -0.40;
    h = 0.70;
    
    %% Perdidas de la primera y segunda pareja de engranajes
    
    K = (1000*Tm*(z1+z2)) / (2*bw*(rw_1)^2*z2); %Factor de carga
    fm = (v^j * K^g) / (C1 * V_pitch^h);        %Coeficiente de fricción
    
    Hs = (u+1)*(( ((ro_2)^2 / (rw_2)^2) - cosd(alpha)^2 )^0.5 - sind(alpha)); %Relacion de deslizamiento al principio de la aproximación 
    Ht = ((u+1)/u)*(( ((ro_1)^2 / (rw_1)^2) - cosd(alpha)^2 )^0.5 - sind(alpha)); %Relacion de deslizamiento al final del receso
    
    M = (2*cosd(alpha)*(Hs + Ht)) / (Hs^2 + Ht^2); %Ventaja mecánica del engranado
    
    Pm = (fm*Tm*n1*cosd(betha)^2) / (9549*M); %Perdidas de potencia por el engranado en kW

end


