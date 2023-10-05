%% Restricciones optimizacion 
function [g_ineq, h_eq]= moo_const (x)
    
    %Se renombran las variables a iterar
    I_primera=x(1); I_segunda=x(2);
    I_E=3.73; %Constante relación transmisión final-optimizar torque (Ver General)

    global HB St S_ut Rpm_max Tor_max z11 z12 z21 z22 m d1e d2e R_rodadura R_pendiente R_inercia R_aerodinamica rueda_carga Fr1
       
    [z11, z12, z21, z22, m, U1, Fr1, Fa1, U2, Fr2, Fa2]=Calculo_engranajes(I_primera, I_segunda, HB, Rpm_max, Tor_max, St);
    [z1d, z2d, m_d, U_d, Fr_d, Fa_d]=Calculo_diferencial(I_E, Tor_max, Rpm_max, I_primera, St);
    [d1e]=Eje_principal(U1, Fr1, Fa1, z11, U2, Fr2, Fa2, z21, m, S_ut);
    [d2e]=Eje_secundario(U1, Fr1, Fa1, z12, U2, Fr2, Fa2, z22, m, S_ut);
    
    % Se obtiene el torque requerido para mover el vehiculo 
    Coeficiente_seguridad=0.95; 
    
    R_resistente1=R_inercia+R_rodadura+R_pendiente;
    R_resistente2=R_rodadura+R_aerodinamica;

    T_resistente1=R_resistente1*rueda_carga;
    T_resistente2=R_resistente2*rueda_carga;
    
    % Se obtiene el torque generador por la transmisión en las ruedas
    T_rueda_primera=(Tor_max*Coeficiente_seguridad)/(1/(I_E*I_primera));
    T_rueda_segunda=(Tor_max*Coeficiente_seguridad)/(1/(I_E*I_segunda));
    
    %% Restricciones de desigualdad
    
    % Esfuerzo admisible mayor al esfuerzo en engranajes
    % Número de dientes mínimo en los engranajes (evitar interferencias)
    % Torque resistente menor al torque generado en la 2° relación
    % Torque resistente menor al torque generado en la 1° relación

    %g_ineq(1)=sigma_f1-(sigma_adm); g_ineq(2)=sigma_f2-sigma_adm;
    g_ineq(1)=16-z12; g_ineq(2)=16-z22; g_ineq(3)=16-z21;
    g_ineq(4)=T_resistente2-T_rueda_segunda;
    g_ineq(5)=T_resistente1-T_rueda_primera;
    
    %% No se plantean restricciones de igualdad
    h_eq=[];

end

%% Funciones utilizadas para obtener las restricciones

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

function [z1d, z2d, m_d, U_d, Fr_d, Fa_d]=Calculo_diferencial(I_E, Tor_max, Rpm_max, I_primera, St)
    
    modulo = 3; %valor asumido estandar de tabla para paso fino en mm
    pd = 1/modulo; %[mm^-1] paso diametral
    %pd = 25.4/modulo; %pulg^-1 paso diametral

    psi = 23; %angulo de helice
    alpha = 20; %angulo de presion   

    T1 = Tor_max*I_primera; % torque del piñon conductor en Nm (en el eje 2, asumiendo primera marcha)
    T_p = T1; %torque aplicado al piñon

    n1 = Rpm_max/I_primera; % rpm del piñon conductor (en el eje 2)

    %% Datos de los engranajes
    z1d = 25;  %Se asume un número de dientes del piñon
    z2d = ceil(I_E*z1d); %Número de dientes del engrane

    modulo_c = modulo/cosd(psi);  %Modulo aparente
    dp_1 = z1d*modulo_c; %Diametro de paso piñon mm
    dp_2 = z2d*modulo_c; %Diameto de paso engrane mm
    b = 12/pd;           %Ancho de la cara en mm
    h_c = 1/pd;            %Addendum
    h_b = 1.25/pd;         %Deddendum
    h = h_c+h_b;          %Profundidad del diente

    pc_1 = (pi*dp_1)/z1d; %Paso circular piñon
    pc_2 = (pi*dp_2)/z2d; %Paso cicular engranaje

    p_t = pi/pd;          %Paso transversal
    p_n = p_t*cosd(psi);  %Paso normal
    p_x = p_n/sind(psi);  %Paso axial
    pd_nd = pd/cosd(psi); %Paso diametral en el plano normal
    alpha_n= atand(cosd(alpha)*tand(psi)); %Angúlo de presion normal 

    R_1 = dp_1/2;    R_2 = dp_2/2;
    C = R_1 + R_2;
    a_p = h_c;    a_g = a_p;

    %Longitud de accion en mm
    Z = ( (R_1 + a_p )^2 - (R_1*cosd(alpha))^2 )^0.5 + ( (R_2 + a_g )^2 - (R_2*cosd(alpha))^2 )^0.5 - C*sind(alpha); 
    mp = (pd*Z)/(pi*cosd(alpha));  %Razón de contacto transversal
    mf = (b*pd*tand(psi))/pi;   %Razón de contacto axial, debe ser mayor a 1.15

    %% Fuerzas en el diferencial

    U_d = (T_p*1000)/R_1;  %fuerza tangencial en N
    Fr_d = U_d*tand(alpha); %fuerza radial en N
    Fa_d = U_d*tand(psi); %fuerza axial en N
    %W_total = U_d/(cosd(psi)*cosd(alpha_n)); %fuerza total en N

    %% Factores de comportamiento
    
    %Factor dinamico
    Qv = 11; %Indice de calidad del engranaje
    B = ((12-Qv)^(2/3))/4;
    A = 50 + 56*(1-B);

    Vt = R_1*((2*pi)/(60*1000))*n1;
    Kv = (A/(A + (200*Vt)^0.5 ) )^B;

    % Factor de distribucion de carga
    Km = 1.6;

    % Factor aplicacion, fuente de energia uniforme y carga con impacto moderado
    Ka = 1.25;

    % Factor de tamaño según el modulo
    Ks=2E-05*modulo^4 - 0.0011*modulo^3 + 0.0164*modulo^2 - 0.0647*modulo + 1.06;

    Kl=1;  %Factor engranaje loco
    Kb=1;  %Factor espesor aro para engranaje solido

    % Factor de Lewis
    Y = -0.00000005*z1d^4+0.00001*z1d^3-0.0007*z1d^2+0.0241*z1d+0.05;

    %% Esfuerzos
    sigma_adm=0;    sigma_b=1;
    NS=1.5; %Factor de seguridad

    while (sigma_adm/sigma_b<NS)
    sigma_b = (U_d/(b*modulo*Y))*((Ka*Km*Ks*Kb*Kl)/Kv); %esfuerzo de flexion en el diente en MPa

    % Asumiendo temperatura en caja de 80° C
    kt=(460+(80*9/5)+32)/620;   %Factor temperatura
    kr=1;                    %Factor fiabilidad    
    kn=9.4518*(10^6)^(-0.148);  %Factor duracion
    sigma_adm=St*kn/(kt*kr); %Esfuerzo admisible   

        if (sigma_adm/sigma_b<NS)
            modulo=modulo+0.5;  
        end
    end

    m_d=modulo;

end

function [d1e]=Eje_principal(U1, Fr1, Fa1, z11, U2, Fr2, Fa2, z21, m, S_ut)
    
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
    if (roundn(d_1m,-1)~=roundn(d,-1))
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

function [d2e]=Eje_secundario(U1, Fr1, Fa1, z12, U2, Fr2, Fa2, z22, m, S_ut)
 
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
