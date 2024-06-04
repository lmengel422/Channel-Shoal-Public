% ***************************************************************************
%  Time advancement code
%   Steps a single timestep for u, v, c, rho, q2, Q2L2, l, kz, nu_t, kq
%   All diffusion/viscous terms handled implicitly
% ***************************************************************************

%% FIRST WATER COLUMN--------------------------------------------------------
if first_water_column == 1
    % ***************************************************************************
    %  Update pressure forcing term for the current timestep
    % ***************************************************************************
    for i = 1:N
        if T_Px == 0.0
            Px(i) = Px0;  %#ok<*SAGROW> % Steady and constant forcing for now
        else
            Px(i) = Px0*(1+springneap*cos(2*pi*t(m)/(3600.*24*14)))*cos(2*pi*t(m)/(3600.*T_Px))+g*alpha*Cx*z(i)+coeff*g*alpha*Cx*H; %Added baroclinic, freshwater creep, springneap
        end
    end

    %Update shear velocity at bottom boundary for use later
    ustar = abs(U(1))*sqrt(C_D); % Explicit dependence on C_D

    % ****************************************************************************
    % Update parameters for the model, Sm and Sh
    % ****************************************************************************
    for i=1:N		
       Gh=-(N_BV(i)*L(i)/(Q(i)+SMALL))^2; 
       %set LIMITER for Gh 
       Gh=min(Gh, 0.0233);
       Gh=max(Gh, -0.28);
       num=B1^(-1/3)-A1*A2*Gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1));
       dem=(1-3*A2*Gh*(B2+6*A1))*(1-9*A1*A2*Gh);
       Sm(i)=num/dem;
       Sh(i)=A2*(1-6*A1/B1)/(1-3*A2*Gh*(B2+6*A1)); 
    end

    % ***************************************************************************
    %  Place previous variable f into fp (i.e. u into up, q2 into q2p, etc)
    % ***************************************************************************
    Up = U;
    Vp = V;
    Cp = C;
    Biop = Bio;
    Q2p = Q2;
    Q2Lp = Q2L;
    Lp = L;
    Kzp = Kz;
    nu_tp = nu_t;
    Kqp = Kq;
    N_BVp = N_BV;

    % ***************************************************************************
    %  Advance velocity (U,V)
    % ***************************************************************************
    for i = 2:N-1
       aU(i) = -0.5*beta*(nu_tp(i)+nu_tp(i-1));
       bU(i) = 1+0.5*beta*(nu_tp(i+1)+2*nu_tp(i)+nu_tp(i-1));
       cU(i) = -0.5*beta*(nu_tp(i)+nu_tp(i+1));
       dU(i) = Up(i) - dt*Px(i);
    end
    %bottom boundary: log-law
    bU(1) = 1+0.5*beta*(nu_tp(2)+nu_tp(1)+2*(sqrt(C_D)/kappa)*nu_tp(1));
    cU(1) = -0.5*beta*(nu_tp(2)+nu_tp(1));
    dU(1) = Up(1) - dt*Px(1);
    %top boundary: no stress
    aU(N) = -0.5*beta*(nu_tp(N)+nu_tp(N-1));
    bU(N) = 1+0.5*beta*(nu_tp(N)+nu_tp(N-1));
    dU(N) = Up(N) - dt*Px(N);
    % Use Thomas algorithm to solve for U
    for i=2:N
       bU(i) = bU(i) - aU(i)/bU(i-1)*cU(i-1);
       dU(i) = dU(i) - aU(i)/bU(i-1)*dU(i-1);
    end
    U(N) = dU(N)/bU(N);
    for i=N-1:-1:1
      U(i) = 1/bU(i)*(dU(i)-cU(i)*U(i+1));
    end

    % ***************************************************************************
    %  Advance scalars/density (C, rho, Bio) 
    % ***************************************************************************
    %Salinity 
    for i = 2:N-1
       aC(i) = -0.5*beta*(Kzp(i)+Kzp(i-1));
       bC(i) = 1+0.5*beta*(Kzp(i+1)+2*Kzp(i)+Kzp(i-1));
       cC(i) = -0.5*beta*(Kzp(i)+Kzp(i+1));
       if i>N-N2
           dC(i) = Cp(i)-dt*U(i)*Cx+dt*ratio1*exchange(2,i-(N-N2)); 
       else
           dC(i) = Cp(i)-dt*U(i)*Cx;
       end
    end
    %bottom boundary: no flux
    bC(1) = 1+0.5*beta*(Kzp(2)+Kzp(1));
    cC(1) = -0.5*beta*(Kzp(2)+Kzp(1));
    dC(1) = Cp(1)-dt*U(1)*Cx;
    %top boundary: no flux for scalars
    aC(N) = -0.5*beta*(Kzp(N)+Kzp(N-1));
    bC(N) = 1+0.5*beta*(Kzp(N)+Kzp(N-1));
    dC(N) = Cp(N)-dt*U(N)*Cx+dt*ratio1*exchange(2,N2); 
    % Use Thomas algorithm to solve for C
    for i=2:N
       bC(i) = bC(i) - aC(i)/bC(i-1)*cC(i-1);
       dC(i) = dC(i) - aC(i)/bC(i-1)*dC(i-1);
    end
    C(N) = dC(N)/bC(N);
    C(N) = max(C(N),min_C);    %Boundary conditions
    C(N) = min(C(N),max_C); 
    for i=N-1:-1:1
     C(i) = 1/bC(i)*(dC(i)-cC(i)*C(i+1));
     C(i) = max(C(i),min_C);    %Boundary conditions
     C(i) = min(C(i),max_C); 
    end
    
    %update density and Brunt-Vaisala frequency
    for i = 1:N
       rho(i)=rho0*(1-alpha*(C(i)-T)); % Single scalar, linear equation of state
    end
    N_BV(1)=sqrt((-g/rho0)*(rho(2)-rho(1))/(dz));
    for i =2:N-1
       N_BV(i)=sqrt((-g/rho0)*(rho(i+1)-rho(i-1))/(2*dz));
    end
    N_BV(N)=sqrt((-g/rho0)*(rho(N)-rho(N-1))/(dz));

    %Update irradiance for option 1
    if irradiance == 1
        if irr_type == 2
            kb = kb_coeff*cumsum(Biop); %[m-1] Self-shading by phyto-plankton biomass (correct way)
        else
            kb = kb_coeff*Biop; %[m-1] Self-shading by phyto-plankton biomass
        end
        if irr_type == 3
            %I = I0.*exp((-kt-kb).*abs(z));
            I(1) = I0*exp(-(kt+kb(1))*dz/2);
            for i=2:N
                I(i)=I(i-1)*exp(-(kt+kb(i))*dz);
            end
            I_flipped = flip(I); %Flip I so high irradiance at the top
        else
            I = I0.*exp((-kt).*abs(z)-kb);
        end
       % f_I = tanh(a*I); 
        f_I = tanh(a*I_flipped); 
        munet = Pmax*(f_I-r)/theta;
    end

    %Biomass 
    for i = 2:N-1
       aBio(i) = -0.5*beta*(Kzp(i)+Kzp(i-1)); 
       bBio(i) = 1+0.5*beta*(Kzp(i+1)+2*Kzp(i)+Kzp(i-1))-dt*(munet(i)-ZP)+dt/dz*ws; %Upwind settling velocity 
       cBio(i) = -0.5*beta*(Kzp(i)+Kzp(i+1))-dt/dz*ws;
       if i>N-N2
           dBio(i) = Biop(i)-dt*U(i)*Biox+dt*ratio1*exchange(1,i-(N-N2)); 
       else
           dBio(i) = Biop(i)-dt*U(i)*Biox;
       end
    end
    %bottom boundary: settling (BG flux)
    bBio(1) = 1+0.5*beta*(Kzp(2)+Kzp(1))-dt*(munet(1)-ZP)+dt/dz*ws;
    cBio(1) = -0.5*beta*(Kzp(2)+Kzp(1))-dt/dz*ws;
    F_R = -(BG*Biop(1)); %Bottom benthic grazing flux
    dBio(1) = Biop(1)-dt*U(1)*Biox+dt/dz*F_R; %No exchange since bottom of channel
    %top boundary: no flux for scalars
    aBio(N) = -0.5*beta*(Kzp(N)+Kzp(N-1));
    bBio(N) = 1+0.5*beta*(Kzp(N)+Kzp(N-1))-dt*(munet(N)-ZP)+dt/dz*ws;
    dBio(N) = Biop(N)-dt*U(N)*Biox+dt*ratio1*exchange(1,N2);
    % Use Thomas algorithm to solve for C
    for i=2:N
       bBio(i) = bBio(i) - aBio(i)/bBio(i-1)*cBio(i-1);
       dBio(i) = dBio(i) - aBio(i)/bBio(i-1)*dBio(i-1);
    end
    Bio(N) = dBio(N)/bBio(N);
    Bio(N) = max(Bio(N),min_Bio);   %Boundary conditions
    for i=N-1:-1:1
      Bio(i) = 1/bBio(i)*(dBio(i)-cBio(i)*Bio(i+1));
      Bio(i) = max(Bio(i),min_Bio);    %Boundary conditions
    end

    % ***************************************************************************
    %  Advance turbulence parameters (q2, q2l - q2 first, then q2l)
    % ***************************************************************************
    for i=2:N-1   
       diss = 2*dt*(Q2p(i))^(1/2)/(B1*Lp(i)); % coefficient for linearized term
       aQ2(i) = -0.5*beta*(Kqp(i)+Kqp(i-1));
       bQ2(i) = 1+0.5*beta*(Kqp(i+1)+2*Kqp(i)+Kqp(i-1))+diss;
       cQ2(i) = -0.5*beta*(Kqp(i)+Kqp(i+1));
       dQ2(i) = Q2p(i) + 0.25*beta*(nu_tp(i))*(Up(i+1)-Up(i-1))^2 - dt*(Kzp(i))*(N_BVp(i))^2;
    end
    %Bottom boundary (i=1)
    Q2bot=B1^(2/3)*ustar^2; %Q2(0) is prescribed
    bndryterm = 0.5*beta*(Kqp(1))*Q2bot;
    diss = 2*dt*(Q2p(1))^(1/2)/(B1*Lp(1));
    bQ2(1) = 1+0.5*beta*(Kqp(2)+Kqp(1))+diss;
    cQ2(1) = -0.5*beta*(Kqp(1)+Kqp(2));
    dQ2(1) = Q2p(1) + dt*(ustar^4)/nu_tp(i)- dt*(Kzp(1))*(N_BVp(1))^2 + bndryterm;
    %Top boundary (i=N)
    diss = 2*dt*(Q2p(N))^(1/2)/(B1*Lp(N));
    aQ2(N) = -0.5*beta*(Kqp(N)+Kqp(N-1));
    bQ2(N) = 1+0.5*beta*(Kqp(N)+2*Kqp(N)+Kq(N-1))+diss;
    dQ2(N) = Q2p(N) + 0.25*beta*(nu_tp(N))*(Up(N)-Up(N-1))^2 - 4*dt*(Kzp(N))*(N_BVp(N))^2 ;
    % Use Thomas algorithm to solve for Q2
    for i=2:N
       bQ2(i) = bQ2(i) - aQ2(i)/bQ2(i-1)*cQ2(i-1);
       dQ2(i) = dQ2(i) - aQ2(i)/bQ2(i-1)*dQ2(i-1);
    end
    Q2(N) = dQ2(N)/bQ2(N);
    for i=N-1:-1:1
        Q2(i) = 1/bQ2(i)*(dQ2(i)-cQ2(i)*Q2(i+1));
    end

    %  Kluge to prevent negative values from causing instabilities
    for i=1:N
       if Q2(i)<SMALL
          Q2(i)=SMALL;
       end
    end
    % *******************************************************************

    for i=2:N-1
       diss=2*dt*(Q2p(i))^(1/2)/(B1*Lp(i))*(1+E2*(Lp(i)/(kappa*abs(-H-z(i))))^2+E3*(Lp(i)/(kappa*abs(z(i))))^2);
       aQ2L(i) = -0.5*beta*(Kqp(i)+Kqp(i-1));
       bQ2L(i) = 1+0.5*beta*(Kqp(i+1)+2*Kqp(i)+Kqp(i-1))+diss;
       cQ2L(i) = -0.5*beta*(Kqp(i)+Kqp(i+1));
       dQ2L(i) = Q2Lp(i) + 0.25*beta*(nu_tp(i))*E1*Lp(i)*(Up(i+1)-Up(i-1))^2 ...
            - 2*dt*Lp(i)*E1*(Kzp(i))*(N_BVp(i))^2;
    end
    %Bottom boundary (i=1)
    Q2Lbot= B1^(2/3)*ustar^2*kappa*zb; 
    bndryterm= 0.5*beta*(Kqp(1))*Q2Lbot;
    diss=2*dt*(Q2p(1))^(1/2)/(B1*Lp(1))*(1+E2*(Lp(1)/(kappa*abs(-H-z(1))))^2+E3*(Lp(1)/(kappa*abs(z(1))))^2);
    bQ2L(1) = 1+0.5*beta*(Kqp(2)+Kqp(1))+diss;
    cQ2L(1) = -0.5*beta*(Kq(1)+Kq(2));
    dQ2L(1) = Q2Lp(1) + dt*((ustar^4)/nu_tp(1))*E1*Lp(1)-dt*Lp(1)*E1*(Kzp(1))*(N_BVp(1))^2+bndryterm;
    %Top boundary (i=N)
    diss=2*dt*(Q2p(N))^(1/2)/(B1*Lp(N))*(1+E2*(Lp(N)/(kappa*abs(-H-z(N))))^2+E3*(Lp(N)/(kappa*abs(z(N))))^2);
    aQ2L(N) = -0.5*beta*(Kqp(N)+Kqp(N-1));
    bQ2L(N) = 1+0.5*beta*(Kqp(N)+2*Kqp(N)+Kqp(N-1))+diss;
    dQ2L(N) = Q2Lp(N) + 0.25*beta*(nu_tp(N))*E1*Lp(N)*(Up(N)-Up(N-1))^2 - 2*dt*Lp(N)*E1*(Kzp(N))*(N_BVp(N))^2;
    % Use Thomas algorithm to solve for Q2L
    for i=2:N
       bQ2L(i) = bQ2L(i) - aQ2L(i)/bQ2L(i-1)*cQ2L(i-1);
       dQ2L(i) = dQ2L(i) - aQ2L(i)/bQ2L(i-1)*dQ2L(i-1);
    end
    Q2L(N) = dQ2L(N)/bQ2L(N);
    for i=N-1:-1:1
       Q2L(i) = 1/bQ2L(i)*(dQ2L(i)-cQ2L(i)*Q2L(i+1));
    end

    %  Kluge to prevent negative values from causing instabilities
    for i=1:N
       if Q2L(i)<SMALL
          Q2L(i)=SMALL;
       end	
    end

    % ***************************************************************************
    %  Calculate turbulent lengthscale (l) and mixing coefficients (kz, nu_t, kq)
    %     Works will all updated values 
    % ***************************************************************************
    for i=1:N		
       Q(i)=sqrt(Q2(i));
       L(i)=Q2L(i)/(Q2(i)+SMALL); 
       %limit due to stable stratification
       if (L(i)^2*(N_BV(i))^2>0.281*Q2(i)) 
            %Adjust Q2L as well as L
            Q2L(i)=Q2(i)*sqrt(0.281*Q2(i)/((N_BV(i))^2+SMALL));
            L(i)=Q2L(i)/Q2(i);
       end
       %Keep L from becoming zero -- zb=bottom roughness parameter
       if  (abs(L(i))<=zb)
            L(i)=zb;
       end
       %update diffusivities
       Kq(i)=Sq*Q(i)*L(i) + nu; 
       nu_t(i)=Sm(i)*Q(i)*L(i) + nu; 
       Kz(i)=Sh(i)*Q(i)*L(i) + nu;
    end
end

%% SECOND WATER COLUMN--------------------------------------------------------
if second_water_column == 1
    % ***************************************************************************
    %  Update pressure forcing term for the current timestep
    % ***************************************************************************
    for i = 1:N2
       if T_Px2 == 0.0
          Px2(i) = Px02;  % Steady and constant forcing for now
       else
          Px2(i) = Px02*(1+springneap2*cos(2*pi*t(m)/(3600.*24*14)))*cos(2*pi*t(m)/(3600.*T_Px2))+g*alpha2*Cx2*z2(i)+coeff2*g*alpha2*Cx2*H2; %Added baroclinic, freshwater creep and springneap
       end
    end

    %Update shear velocity at bottom boundary for use later
    ustar2 = abs(U2(1))*sqrt(C_D); % Explicit dependence on C_D
    % ****************************************************************************
    % Update parameters for the model, Sm and Sh
    % ****************************************************************************

    for i=1:N2		
       Gh2=-(N_BV2(i)*L2(i)/(Q_2(i)+SMALL))^2; 
       %set LIMITER for Gh 
       Gh2=min(Gh2, 0.0233);
       Gh2=max(Gh2, -0.28);
       num2=B12^(-1/3)-A12*A22*Gh2*((B22-3*A22)*(1-6*A12/B12)-3*C12*(B22+6*A12));
       dem2=(1-3*A22*Gh2*(B22+6*A12))*(1-9*A12*A22*Gh2);
       Sm2(i)=num2/dem2;
       Sh2(i)=A22*(1-6*A12/B12)/(1-3*A22*Gh2*(B22+6*A12)); 
    end

    % ***************************************************************************
    %  Place previous variable f into fp (i.e. u into up, q2 into q2p, etc)
    % ***************************************************************************
    Up2 = U2;
    Vp2 = V2;
    Cp2 = C2;
    Biop2 = Bio2;
    Q2p2 = Q22;
    Q2Lp2 = Q2L2;
    Lp2 = L2;
    Kzp2 = Kz2;
    nu_tp2 = nu_t2;
    Kqp2 = Kq2;
    N_BVp2 = N_BV2;

    % ***************************************************************************
    %  Advance velocity (U,V)
    % ***************************************************************************
    for i = 2:N2-1
       aU2(i) = -0.5*beta*(nu_tp2(i)+nu_tp2(i-1));
       bU2(i) = 1+0.5*beta*(nu_tp2(i+1)+2*nu_tp2(i)+nu_tp2(i-1));
       cU2(i) = -0.5*beta*(nu_tp2(i)+nu_tp2(i+1));
       dU2(i) = Up2(i) - dt*Px2(i);
    end
    %bottom boundary: log-law
    bU2(1) = 1+0.5*beta*(nu_tp2(2)+nu_tp2(1)+2*(sqrt(C_D)/kappa)*nu_tp2(1));
    cU2(1) = -0.5*beta*(nu_tp2(2)+nu_tp2(1));
    dU2(1) = Up2(1) - dt*Px2(1);
    %top boundary: no stress
    aU2(N2) = -0.5*beta*(nu_tp2(N2)+nu_tp2(N2-1));
    bU2(N2) = 1+0.5*beta*(nu_tp2(N2)+nu_tp2(N2-1));
    dU2(N2) = Up2(N2) - dt*Px2(N2);
    % Use Thomas algorithm to solve for U
    for i=2:N2
       bU2(i) = bU2(i) - aU2(i)/bU2(i-1)*cU2(i-1);
       dU2(i) = dU2(i) - aU2(i)/bU2(i-1)*dU2(i-1);
    end
    U2(N2) = dU2(N2)/bU2(N2);
    for i=N2-1:-1:1
      U2(i) = 1/bU2(i)*(dU2(i)-cU2(i)*U2(i+1));
    end

    % ***************************************************************************
    %  Advance scalars/density (C, rho, Bio) 
    % ***************************************************************************
    %Salinity 
    for i = 2:N2-1
       aC2(i) = -0.5*beta*(Kzp2(i)+Kzp2(i-1));
       bC2(i) = 1+0.5*beta*(Kzp2(i+1)+2*Kzp2(i)+Kzp2(i-1));
       cC2(i) = -0.5*beta*(Kzp2(i)+Kzp2(i+1));
       dC2(i) = Cp2(i)-dt*U2(i)*Cx2-dt*ratio2*exchange(2,i); 
    end
    %bottom boundary: no flux
    bC2(1) = 1+0.5*beta*(Kzp2(2)+Kzp2(1));
    cC2(1) = -0.5*beta*(Kzp2(2)+Kzp2(1));
    dC2(1) = Cp2(1)-dt*U2(1)*Cx2-dt*ratio2*exchange(2,1); 
    %top boundary: no flux for scalars
    aC2(N2) = -0.5*beta*(Kzp2(N2)+Kzp2(N2-1));
    bC2(N2) = 1+0.5*beta*(Kzp2(N2)+Kzp2(N2-1));
    dC2(N2) = Cp2(N2)-dt*U2(N2)*Cx2-dt*ratio2*exchange(2,N2); 
    % Use Thomas algorithm to solve for C
    for i=2:N2
       bC2(i) = bC2(i) - aC2(i)/bC2(i-1)*cC2(i-1);
       dC2(i) = dC2(i) - aC2(i)/bC2(i-1)*dC2(i-1);
    end
    C2(N2) = dC2(N2)/bC2(N2);
    C2(N2) = max(C2(N2),min_C);    %Boundary conditions
    C2(N2) = min(C2(N2),max_C); 
    for i=N2-1:-1:1
     C2(i) = 1/bC2(i)*(dC2(i)-cC2(i)*C2(i+1));
     C2(i) = max(C2(i),min_C);    %Boundary conditions
     C2(i) = min(C2(i),max_C); 
    end
    
    %update density and Brunt-Vaisala frequency
    for i = 1:N2
       rho2(i)=rho02*(1-alpha2*(C2(i)-T2)); % Single scalar, linear equation of state
    end
    N_BV2(1)=sqrt((-g/rho02)*(rho2(2)-rho2(1))/(dz));
    for i =2:N2-1
       N_BV2(i)=sqrt((-g/rho02)*(rho2(i+1)-rho2(i-1))/(2*dz));
    end
    N_BV2(N2)=sqrt((-g/rho02)*(rho2(N2)-rho2(N2-1))/(dz));
    
    %Update irradiance for option 1
    if irradiance == 1
        if irr_type == 2
            kb2 = kb2_coeff*cumsum(Biop2); %[m-1] Self-shading by phyto-plankton biomass (correct way)
        else
            kb2 = kb2_coeff*Biop2; %[m-1] Self-shading by phyto-plankton biomass
        end
        if irr_type == 3
            %I2 = I02.*exp((-kt2-kb2).*abs(z2));
            I2(1) = I02*exp(-(kt2+kb2(1))*dz/2);
            for i=2:N2
                I2(i)=I2(i-1)*exp(-(kt2+kb2(i))*dz);
            end
            I2_flipped = flip(I2); %Need to flip direction I so that highest at surface, lowest at bottom
        else
            I2 = I02.*exp((-kt2).*abs(z2)-kb2);
        end
       % f_I2 = tanh(a2*I2);
        f_I2 = tanh(a2*I2_flipped);
        munet2 = Pmax2*(f_I2-r2)/theta2; 
    end

    %Biomass
    for i = 2:N2-1
       aBio2(i) = -0.5*beta*(Kzp2(i)+Kzp2(i-1)); 
       bBio2(i) = 1+0.5*beta*(Kzp2(i+1)+2*Kzp2(i)+Kzp2(i-1))-dt*(munet2(i)-ZP2)+dt/dz*ws2; %Upwind settling velocity 
       cBio2(i) = -0.5*beta*(Kzp2(i)+Kzp2(i+1))-dt/dz*ws2;
       dBio2(i) = Biop2(i)-dt*U2(i)*Biox2-dt*ratio2*exchange(1,i); 
    end
    %bottom boundary: settling (BG flux)
    bBio2(1) = 1+0.5*beta*(Kzp2(2)+Kzp2(1))-dt*(munet2(1)-ZP2)+dt/dz*ws2;
    cBio2(1) = -0.5*beta*(Kzp2(2)+Kzp2(1))-dt/dz*ws2;
    F_R2 = -(BG2*Biop2(1)); %Bottom benthic grazing flux
    dBio2(1) = Biop2(1)-dt*U2(1)*Biox2-dt*ratio2*exchange(1,1)+dt/dz*F_R2; %Exchange bottom of shoal, but not channel
    %top boundary: no flux for scalars
    aBio2(N2) = -0.5*beta*(Kzp2(N2)+Kzp2(N2-1));
    bBio2(N2) = 1+0.5*beta*(Kzp2(N2)+Kzp2(N2-1))-dt*(munet2(N2)-ZP2)+dt/dz*ws2; 
    dBio2(N2) = Biop2(N2)-dt*U2(N2)*Biox2-dt*ratio2*exchange(1,N2); 
    % Use Thomas algorithm to solve for C
    for i=2:N2
       bBio2(i) = bBio2(i) - aBio2(i)/bBio2(i-1)*cBio2(i-1);
       dBio2(i) = dBio2(i) - aBio2(i)/bBio2(i-1)*dBio2(i-1);
    end
    Bio2(N2) = dBio2(N2)/bBio2(N2);
    Bio2(N2) = max(Bio2(N2),min_Bio);   %Boundary conditions
    for i=N2-1:-1:1
      Bio2(i) = 1/bBio2(i)*(dBio2(i)-cBio2(i)*Bio2(i+1));
      Bio2(i) = max(Bio2(i),min_Bio);    %Boundary conditions
    end

    % ***************************************************************************
    %  Advance turbulence parameters (q2, q2l - q2 first, then q2l)
    % ***************************************************************************
    for i=2:N2-1   
       diss2 = 2*dt*(Q2p2(i))^(1/2)/(B12*Lp2(i)); % coefficient for linearized term
       aQ22(i) = -0.5*beta*(Kqp2(i)+Kqp2(i-1));
       bQ22(i) = 1+0.5*beta*(Kqp2(i+1)+2*Kqp2(i)+Kqp2(i-1))+diss2;
       cQ22(i) = -0.5*beta*(Kqp2(i)+Kqp2(i+1));
       dQ22(i) = Q2p2(i) + 0.25*beta*(nu_tp2(i))*(Up2(i+1)-Up2(i-1))^2 - dt*(Kzp2(i))*(N_BVp2(i))^2;
    end
    %Bottom boundary (i=1)
    Q2bot2=B12^(2/3)*ustar2^2; %Q2(0) is prescribed
    bndryterm2 = 0.5*beta*(Kqp2(1))*Q2bot2;
    diss2 = 2*dt*(Q2p2(1))^(1/2)/(B12*Lp2(1));
    bQ22(1) = 1+0.5*beta*(Kqp2(2)+Kqp2(1))+diss2;
    cQ22(1) = -0.5*beta*(Kqp2(1)+Kqp2(2));
    dQ22(1) = Q2p2(1) + dt*(ustar2^4)/nu_tp2(i)- dt*(Kzp2(1))*(N_BVp2(1))^2 + bndryterm2;
    %Top boundary (i=N)
    diss2 = 2*dt*(Q2p2(N2))^(1/2)/(B12*Lp2(N2));
    aQ22(N2) = -0.5*beta*(Kqp2(N2)+Kqp2(N2-1));
    bQ22(N2) = 1+0.5*beta*(Kqp2(N2)+2*Kqp2(N2)+Kq2(N2-1))+diss2;
    dQ22(N2) = Q2p2(N2) + 0.25*beta*(nu_tp2(N2))*(Up2(N2)-Up2(N2-1))^2 - 4*dt*(Kzp2(N2))*(N_BVp2(N2))^2 ;
    % Use Thomas algorithm to solve for Q2
    for i=2:N2
       bQ22(i) = bQ22(i) - aQ22(i)/bQ22(i-1)*cQ22(i-1);
       dQ22(i) = dQ22(i) - aQ22(i)/bQ22(i-1)*dQ22(i-1);
    end
    Q22(N2) = dQ22(N2)/bQ22(N2);
    for i=N2-1:-1:1
        Q22(i) = 1/bQ22(i)*(dQ22(i)-cQ22(i)*Q22(i+1));
    end

    %  Kluge to prevent negative values from causing instabilities
    for i=1:N2
       if Q22(i)<SMALL
          Q22(i)=SMALL;
       end
    end
    % *******************************************************************

    for i=2:N2-1
       diss2=2*dt*(Q2p2(i))^(1/2)/(B12*Lp2(i))*(1+E22*(Lp2(i)/(kappa*abs(-H2-z2(i))))^2+E32*(Lp2(i)/(kappa*abs(z2(i))))^2);
       aQ2L2(i) = -0.5*beta*(Kqp2(i)+Kqp2(i-1));
       bQ2L2(i) = 1+0.5*beta*(Kqp2(i+1)+2*Kqp2(i)+Kqp2(i-1))+diss2;
       cQ2L2(i) = -0.5*beta*(Kqp2(i)+Kqp2(i+1));
       dQ2L2(i) = Q2Lp2(i) + 0.25*beta*(nu_tp2(i))*E12*Lp2(i)*(Up2(i+1)-Up2(i-1))^2 ...
            - 2*dt*Lp2(i)*E12*(Kzp2(i))*(N_BVp2(i))^2;
    end
    %Bottom boundary (i=1)
    Q2Lbot2= B12^(2/3)*ustar2^2*kappa*zb; 
    bndryterm2= 0.5*beta*(Kqp2(1))*Q2Lbot2;
    diss2=2*dt*(Q2p2(1))^(1/2)/(B12*Lp2(1))*(1+E22*(Lp2(1)/(kappa*abs(-H2-z2(1))))^2+E32*(Lp2(1)/(kappa*abs(z2(1))))^2);
    bQ2L2(1) = 1+0.5*beta*(Kqp2(2)+Kqp2(1))+diss2;
    cQ2L2(1) = -0.5*beta*(Kq2(1)+Kq2(2));
    dQ2L2(1) = Q2Lp2(1) + dt*((ustar2^4)/nu_tp2(1))*E12*Lp2(1)-dt*Lp2(1)*E12*(Kzp2(1))*(N_BVp2(1))^2+bndryterm2;
    %Top boundary (i=N)
    diss2=2*dt*(Q2p2(N2))^(1/2)/(B12*Lp2(N2))*(1+E22*(Lp2(N2)/(kappa*abs(-H2-z2(N2))))^2+E32*(Lp2(N2)/(kappa*abs(z2(N2))))^2);
    aQ2L2(N2) = -0.5*beta*(Kqp2(N2)+Kqp2(N2-1));
    bQ2L2(N2) = 1+0.5*beta*(Kqp2(N2)+2*Kqp2(N2)+Kqp2(N2-1))+diss2;
    dQ2L2(N2) = Q2Lp2(N2) + 0.25*beta*(nu_tp2(N2))*E12*Lp2(N2)*(Up2(N2)-Up2(N2-1))^2 - 2*dt*Lp2(N2)*E12*(Kzp2(N2))*(N_BVp2(N2))^2;
    % Use Thomas algorithm to solve for Q2L2
    for i=2:N2
       bQ2L2(i) = bQ2L2(i) - aQ2L2(i)/bQ2L2(i-1)*cQ2L2(i-1);
       dQ2L2(i) = dQ2L2(i) - aQ2L2(i)/bQ2L2(i-1)*dQ2L2(i-1);
    end
    Q2L2(N2) = dQ2L2(N2)/bQ2L2(N2);
    for i=N2-1:-1:1
       Q2L2(i) = 1/bQ2L2(i)*(dQ2L2(i)-cQ2L2(i)*Q2L2(i+1));
    end

    %  Kluge to prevent negative values from causing instabilities
    for i=1:N2
       if Q2L2(i)<SMALL
          Q2L2(i)=SMALL;
       end	
    end

    % ***************************************************************************
    %  Calculate turbulent lengthscale (l) and mixing coefficients (kz, nu_t, kq)
    %     works with all updated values 
    % ***************************************************************************
    for i=1:N2		
       Q_2(i)=sqrt(Q22(i));
       L2(i)=Q2L2(i)/(Q22(i)+SMALL); 
       %limit due to stable stratification
       if (L2(i)^2*(N_BV2(i))^2>0.281*Q22(i)) 
            %Adjust Q2L as well as L
            Q2L2(i)=Q22(i)*sqrt(0.281*Q22(i)/((N_BV2(i))^2+SMALL));
            L2(i)=Q2L2(i)/Q22(i);
       end
       %Keep L from becoming zero -- zb=bottom roughness parameter
       if  (abs(L2(i))<=zb)
            L2(i)=zb;
       end
       %update diffusivities
       Kq2(i)=Sq2*Q_2(i)*L2(i) + nu; 
       nu_t2(i)=Sm2(i)*Q_2(i)*L2(i) + nu; 
       Kz2(i)=Sh2(i)*Q_2(i)*L2(i) + nu;
    end
end

% *******************************************************************
%%  Calculate difference and perform exchange
% *******************************************************************
if first_water_column == 1 && second_water_column == 1 %Both water columns made 
    if imposeKy == 0 %Use time varying Ky unless imposing
        if condition(t(m),U)
            Ky = Kyp;
            KyS = KySp;
        else
            Ky = ones(1,N2)*Kymin;
            KyS = ones(1,N2)*Kymin;
        end
    end
    diffU = U2 - U(N-N2+1:N);
    diffC = C2 - C(N-N2+1:N);
    diffBio = Bio2 - Bio(N-N2+1:N);
    diffQ2 = Q22 - Q2(N-N2+1:N);
    diffQ2L = Q2L2 - Q2L(N-N2+1:N);
    diffrho = rho2 - rho(N-N2+1:N);
    diffL = L2 - L(N-N2+1:N);
    diffnu_t = nu_t2 - nu_t(N-N2+1:N);
    diffKz = Kz2 - Kz(N-N2+1:N);
    diffKq = Kq2 - Kq(N-N2+1:N);
    diffN_BV = N_BV2 - N_BV(N-N2+1:N);
    
    %Exchange 
    Y = [diffBio;diffC];
    for j=1:2  %Bio,C
        if j < 2
            exchange(j,:) = Ky.*Y(j,:); %Bio exchange
        else
            exchange(j,:) = KyS.*Y(j,:); %Salinity exchange
        end
    end
end