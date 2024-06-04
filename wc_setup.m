% ***************************************************************************
%  Initialize all profiles and closure parameters 
%      Call once before time loop
%      Sets all forcing: pressure gradients, stresses, etc.
%      Should be used to adjust initial temperature/salinity profiles
%      Velocity initialized to zero
%      Turbulence quantities initialized to "SMALL"; Lengthscale parabolic
% ***************************************************************************
%% Flags for certain parameters/plots
filename = 'Sink High Px Low Kts Flipped Full';  %Which test round we are on 'Extra Runs Ky 48 out phase';
savefile = [pwd,'/',filename]; %Where to save created files
export = 0; %Set to 1 to export filename to excel
allvars = 0; %Set to 1 to save everything in wc_setup file
saveme = 0; %Set to 1 to save to file
savefigme = 0; %Whether or not to save figures. Also works if saveme selected
movie = 0; %Whether or not to make a movie
wantpng = 1; %Save fig as png
wantfig = 0; %Save fig as fig
colorscheme = 'jet'; %Color scheme for 3D plots parula jet pink copper  
wantsamecolor = 1; %Set to 1 if want colorbars synced Bio plots
if saveme == 1
    savewhat = {'Cx','Cx2','H','H2','Kybottom','Kytop','Kytype','KySbottom',...
        'KyStop','KyStype','Bio0','Bio02','BG','BG2','Biox','Biox2','T','T2',...
        'T_Px','T_Px2','coeff','coeff2','do_coeff','want_self_shade',...
        'pass_scalar','Px0','Px02','ratio1','ratio2','steady','unstratified',...
        'ws','ws2','condnum','kt','kt2','Umax','numdays','ebbfl','Kymin','irr_type'}; 
end
print = 0; %Set to 1 to print parameters used and Si
si_type = 1; %Set to 1 for method 1 Simpson Number, anything else option 2
calc_rous = 0; %Set to 1 to calculate the Rouse number
%% Which plots wanted 
plotme(1) = 0; %Set to 1 to plot WC 1 Salinity
plotme(2) = 0; %Set to 1 to plot WC 1 Velocity
plotme(3) = 0; %Set to 1 to plot WC 1 Turbulent Kinetic Energy
plotme(4) = 0; %Set to 1 to plot WC 1 Bio surfplot
plotme(5) = 0; %Set to 1 to plot WC 1 Bio subplots
plotme(6) = 0; %Set to 1 to plot WC 2 Salinity
plotme(7) = 0; %Set to 1 to plot WC 2 Velocity
plotme(8) = 0; %Set to 1 to plot WC 2 Turbulent Kinetic Energy
plotme(9) = 0; %Set to 1 to plot WC 2 Bio surfplot
plotme(10) = 0; %Set to 1 to plot WC 2 Bio subplots
plotme(11) = 0; %Set to 1 to plot Both WCs Salinity
plotme(12) = 0; %Set to 1 to plot Both WCs Velocity
plotme(13) = 0; %Set to 1 to plot Both WCs TKE
plotme(14) = 0; %Set to 1 to plot Both WCs Bio surfplots
plotme(15) = 0; %Set to 1 to plot Both WCs Bio subplots
plotme(16) = 0; %Set to 1 to plot Both WCs Middle water column
plotme(17) = 1; %Set to 1 to plot Both WCs Depth Averaged Bio
plotme(18) = 0; %Set to 1 to plot separate fluxes of Bio expression
plotme(19) = 0; %Set to 1 to plot WC1 Bottom Boundary
plotme(20) = 0; %Set to 1 to plot WC2 Bottom Boundary
plotme(21) = 0; %Set to 1 to plot WC 1 Eddy Diffusivity
plotme(22) = 0; %Set to 1 to plot WC 2 Eddy Diffusivity
plotme(23) = 0; %Set to 1 to plot Ky square waves
plotme(24) = 0; %Set to 1 to plot scatter plot of binary results against BG, BG2 
plotme(25) = 0; %Set to 1 to plot scatter plot binary results against Kytest,condnum
plotme(26) = 0; %Set to 1 to plot bar graph for each math calculated term WC1
plotme(27) = 0; %Set to 1 to plot bar graph for each math calculated term WC2
plotme(28) = 0; %Set to 1 to plot scatter plot each math calculated term WC1
plotme(29) = 0; %Set to 1 to plot scatter plot each math calculated term WC2
plotme(30) = 0; %Set to 1 to plot model line WC1
plotme(31) = 0; %Set to 1 to plot bar graphs relative
plotme(32) = 0; %Set to 1 to plot contour plots against Kymax, Numdays
plotme(33) = 0; %Set to 1 to plot both WCS depth averaged salinity
plotme(34) = 0; %Set to 1 to plot both WCs depth averaged velocity
plotme(35) = 0; %Set to 2 to rerun code with different Ky and impose onto previous plot. Not set to 1 so doesn't save blanks (only save every other one)
plotme(36) = 0; %Set to 1 to plot surface biomass difference
% ***************************************************************************
%%  Physical parameters for both water columns
% ***************************************************************************
steady = 0; %Set to 1 for steady (T_Px = 0)
unstratified = 1; %Set to 1 to have unstratified (delC = 0)
irradiance = 1; %Set to 1 to perform Lucas irradiance 2 to perform Roy irradiance 
irr_type = 3; %Set to 1 for way done in original runs, set to 2 for cumsum, set to 3 for Lucas way
pass_scalar = 0; %Set to 1 for passive scalar case (alpha=0)
do_coeff = 1; %Set to 1 to have freshwater creep
want_self_shade = 1; %Set to 1 to enable self-shading in both watercolumns I exp

z0=0.01; %bottom roughness [m]
zb=10*z0; %bottom height
g=9.81; %m^2/s - gravity
C_D = 0.0025; %friction coefficient
SMALL=1e-6;
kappa=0.4; %von Karman constant
nu=1e-6; %m^2/s kinematic viscosity

%From Bio Model
min_C=SMALL; %[psu] 
max_C=35; %[psu]
min_Bio=SMALL; %[mg chl a/m^3]
Ky_range = [0,6,7.5,8,10,35,40,48,50]; %[m^2/s] Values Ky to test. Indexed by Kytest
Kymin_range = [0,4,2.5,2,0,15,10,2,0]; %[m^2/s) Bottom value in time varying Ky cycle 
tide_range = [9.1667e-05,0.00015]; %[m/s^2] Linear range from neap to spring tide. 
springneap = 0.25; %Percent influence spring-neap cycle WC1 0.2
springneap2 = 0.25; %Percent influence spring-neap cycle WC2 0.2

%% FIRST WATER COLUMN-------------------------------------------------------
if first_water_column == 1
    % ***************************************************************************
    %  Forcing parameters - Need to modify to allow for time variable Px.
    % ***************************************************************************
    Px0 = mean(tide_range); %[kPa*m^2/1000kg=m/s^2] Magnitude on pressure gradient forcing. 1/rho0 already included.  0.00015 0.00015/2 0.001 0.00010799
    
    if steady == 1
        T_Px = 0; %Set to 0 for steady
    else
        T_Px = 12.4; %Period (hours) on pressure gradient forcing. 
    end
    
    Cx = 0.0005; %[psu/m] Along channel salinity gradient 0.0005 0.0001 
    Biox = 0; %[mg chl a/m^4] Set these to 0 for now. Along channel biology gradient
    
    % ***************************************************************************
    %  Turbulence closure parameters 
    % ***************************************************************************
    A1=0.92;
    A2=0.74;
    B1=16.6;
    B2=10.1;
    C1=0.08;
    E1=1.8;
    E2= 1.33;
    E3=0.25;
    Sq=0.2;
    
    % ********************************************************************
    %  Setup initial conditions for scalar and density
    % ********************************************************************
    if unstratified == 1
        delC=0; %set to zero for Unstratified Case 
        delBio=0; %set to zero for Unstratified Case
    else
        delC=-1; %change in salinity at initial salocline [psu]; Want less saline at top 
        delBio=0; %change in biomass at initial salocline [mg chl a/m^3]; Don't stratify at initial
    end
    
    zdelC = -5; %position of initial salocline 
    dzdelC = -2; %width of initial salocline
    zdelBio = -5; %position of initial salocline
    dzdelBio = -2; %width of initial salocline
    
    if pass_scalar == 1
        alpha = 0; %set to zero for passive scalar case
    else
        alpha = -7.5*10^-4; %coefficient haline contraction [psu-1]
    end
    
    if do_coeff == 1
        coeff = 0.25; %ratio of constant freshwater creep 0.3 
    else
        coeff = 0;
    end
    rho0=1000; %kg/m^3 - water density
    T = 15; %psu scalar salinity. 15 default case
    
    %Salinity set up
    for i=1:N
        C(i)=T; %#ok<*SAGROW> %Initially uniform salinity profile
        if z(i)<=zdelC-0.5*dzdelC
            C(i)=T;
        elseif z(i)>=zdelC+0.5*dzdelC
            C(i)=T+delC;
        else
            C(i)=T+delC*(z(i)-zdelC+0.5*dzdelC)/dzdelC;
        end
       rho(i)=rho0*(1-alpha*(C(i)-T)); % Single scalar, linear equation of state
    end
    
    %Brunt-Vaisala frequency from density profile
    N_BV(1)=sqrt((-g/rho0)*(rho(2)-rho(1))/(dz));
    for i =2:N-1
       N_BV(i)=sqrt((-g/rho0)*(rho(i+1)-rho(i-1))/(2*dz));
    end
    N_BV(N)=sqrt((-g/rho0)*(rho(N)-rho(N-1))/(dz));
    
    %Biomass set up
    Bio0 = 3; %[mg chl a/m^3] Initial biomass 
    for i=1:N
        Bio(i)=Bio0; %Initially uniform salinity profile
        if z(i)<=zdelBio-0.5*dzdelBio
            Bio(i)=Bio0;
        elseif z(i)>=zdelBio+0.5*dzdelBio
            Bio(i)=Bio0+delBio;
        else
            Bio(i)=Bio0+delBio*(z(i)-zdelBio+0.5*dzdelBio)/dzdelBio;
        end
    end
    
    %Biological Parameters
    ws = 0.5/86400; %[m/s] Phytoplankton sinking rate. Negative built in code.
    a = 0.1*86400; %[m^2*s/mol quanta] Photosynthetic efficiency at low irradiance
    I0 = 30/86400; %[mol quanta m^2/s] Mean daily surface irradiance
    Pmax = 100/86400; %[mgC(mg cla)^-1 s^-1] Maximum phytoplankton carbon assimilation rate
    r = 0.05; %[% of Pmax] Respiration rate 
    ZP = 0.1/86400; %[s^-1] Zooplankton grazing rate
    BG = 2/86400; %Fix for now. %[m^3 m^-2 s^-1] Benthic grazing rate 2
    theta = 50; %[mgC(mg chl a)^-1] Ratio phytoplankton cellular carbon to chlorophyll
    kt = 2; %[m-1] Abiotic light attenuation coefficient. Range 1-7
    if want_self_shade == 1
        kb_coeff = 0.016; %Coefficient for self-shading 
    else
        kb_coeff = 0;
    end
    if irr_type == 2
        kb = kb_coeff*cumsum(Bio); %[m-1] Self-shading by phyto-plankton biomass 
    else
        kb = kb_coeff*Bio; %[m-1] Self-shading by phyto-plankton biomass
    end
    
    if irradiance == 1 %If want to perform irradiance term
        if irr_type == 3
            %I = I0.*exp((-kt-kb).*abs(z));
            I(1) = I0*exp(-(kt+kb(1))*dz/2);
            for i=2:N
                I(i)=I(i-1)*exp(-(kt+kb(i))*dz);
            end
            I_flipped = flip(I); %Flip so that high irradiance at top, low at bottom
        else
            I = I0.*exp((-kt).*abs(z)-kb);
        end
     %   f_I = tanh(a*I); 
        f_I = tanh(a*I_flipped); 
    elseif irradiance == 2 %Roy Irradiance
        %Irradiance - parameters from Roy
        Vol = 1.6; %Average body volume phytoplankton
        Ec = 0.2; %Extinction coefficient water
        Ext = 0.12*Vol^-0.33; %Self shading of phytoplankton due to body size
        Iopt = 680; %Optimum surface solar radiation for photosynthesis
        Ir = 800; %Surface solar irradiance.

        %Roy irradiance effect
        I = Ir.*exp((-Ec-Ext)*abs(z)); 
        
        f_I = I./Iopt.*exp(1-I./Iopt);
    else
        f_I = ones(1,N);
    end
    
    munet = Pmax*(f_I-r)/theta; 
    
    % *******************************************************************
    %  Set initial conditions for u, q^2, and other turbulence quantities
    % *******************************************************************
    t(1)=0;
    for i=1:N
       U(i) = 0.0;
       V(i) = 0.0;
       Q2(i)=SMALL;  %"seed" the turbulent field with small values, then let it evolve
       Q2L(i)=SMALL;
       Q(i)=sqrt(Q2(i));
       L(i)=-kappa*H*(z(i)/H)*(1-(z(i)/H)); %Q2L(n,1)/Q2(n,1) = 1 at initialization;
       Gh= -(N_BV(i)*L(i)/(Q(i)+SMALL))^2; 
       Gh=min(Gh,0.0233);
       Gh=max(Gh,-0.28);
       num=B1^(-1/3)-A1*A2*Gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1));
       dem=(1-3*A2*Gh*(B2+6*A1))*(1-9*A1*A2*Gh);
       Sm(i)=num/dem;
       Sh(i)=A2*(1-6*A1/B1)/(1-3*A2*Gh*(B2+6*A1)); 
       Kq(i)=Sq*Q(i)*L(i)+nu;  %turbulent diffusivity for Q2
       nu_t(i)=Sm(i)*Q(i)*L(i)+nu; %turbulent viscosity
       Kz(i)=Sh(i)*Q(i)*L(i)+nu; %turbulent scalar diffusivity
    end
        
    % *******************************************************************
    %  Save initial conditions as first columns in saved matrix
    % *******************************************************************
    Um(:,1) = U; %Velocity
    Cm(:,1) = C; %Salinity
    Biom(:,1) = Bio; %Biomass
    Q2m(:,1) = Q2; %TKE
    Q2Lm(:,1) = Q2L; %TKE*Length scale
    rhom(:,1) = rho; %Density
    Lm(:,1) = L; %Length scale
    nu_tm(:,1) = nu_t; %Eddy viscosity
    Kzm(:,1) = Kz; %Eddy diffusivity
    Kqm(:,1) = Kq; %Effective diffusion coefficient TKE
    N_BVm(:,1) = N_BV; %Brunt-Vaisala frequency
end

%% SECOND WATER COLUMN-------------------------------------------------------
if second_water_column == 1  
    % ***************************************************************************
    %  Forcing parameters - Need to modify to allow for time variable Px.
    % ***************************************************************************
    Px02 = mean(tide_range); %[kPa*m^2/1000kg=m/s^2] Magnitude on pressure gradient forcing. 1/rho0 already included.  0.00015 0.001 0.00010799
    
    if steady == 1
        T_Px2 = 0; %Set to 0 for steady
    else
        T_Px2 = 12.4; %Period (hours) on pressure gradient forcing. 
    end

    Cx2 = 0.0005; %[psu/m] Along channel salinity gradient 0.0005 0.0001 0.00005
    Biox2 = 0; %[mg chl a/m^4] Set these to 0 for now

    % ***************************************************************************
    %  Turbulence closure parameters 
    % ***************************************************************************
    A12=0.92;
    A22=0.74;
    B12=16.6;
    B22=10.1;
    C12=0.08;
    E12=1.8;
    E22= 1.33;
    E32=0.25;
    Sq2=0.2;

    % ********************************************************************
    %  Setup initial conditions for scalar and density
    % ********************************************************************
    if unstratified == 1
        delC2=0; %set to zero for Unstratified Case 
        delBio2=0; %set to zero for Unstratified Case
    else
        delC2=1; %change in salinity at initial salocline [psu]; 
        delBio2=1; %change in biomass at initial salocline [mg chl a/m^3];  
    end
    
    zdelC2 = -5; %position of initial salocline
    dzdelC2 = -2; %width of initial salocline
    zdelBio2 = -5; %position of initial salocline
    dzdelBio2 = -2; %width of initial salocline
    
    if pass_scalar == 1
        alpha2 = 0; %set to zero for passive scalar case
    else
        alpha2 = -7.5*10^-4; %thermal expansivity or coefficient haline contraction [psu-1]
    end
    
    if do_coeff == 1
        coeff2 = 0.3; %ratio of constant freshwater creep 0.3
    else
        coeff2 = 0;
    end
    rho02=1000; %kg/m^3 - water density
    T2 = 15; %psu scalar salinity. 15 default case

    %Salinity set up
    for i=1:N2
        C2(i)=T2; %Initially uniform salinity profile
        if z2(i)<=zdelC2-0.5*dzdelC2
            C2(i)=T2;
        elseif z2(i)>=zdelC2+0.5*dzdelC2
            C2(i)=T2+delC2;
        else
            C2(i)=T2+delC2*(z2(i)-zdelC2+0.5*dzdelC2)/dzdelC2;
        end
       rho2(i)=rho02*(1-alpha2*(C2(i)-T2)); % Single scalar, linear equation of state
    end

    %Brunt-Vaisala frequency from density profile
    N_BV2(1)=sqrt((-g/rho02)*(rho2(2)-rho2(1))/(dz));
    for i =2:N2-1
       N_BV2(i)=sqrt((-g/rho02)*(rho2(i+1)-rho2(i-1))/(2*dz));
    end
    N_BV2(N2)=sqrt((-g/rho02)*(rho2(N2)-rho2(N2-1))/(dz));

    %Biomass set up
    Bio02 = 3; %[mg chl a/m^3] Initial biomass 
    for i=1:N2
        Bio2(i)=Bio02; %Initially uniform salinity profile
        if z2(i)<=zdelBio2-0.5*dzdelBio2
            Bio2(i)=Bio02;
        elseif z2(i)>=zdelBio2+0.5*dzdelBio2
            Bio2(i)=Bio02+delBio2;
        else
            Bio2(i)=Bio02+delBio2*(z2(i)-zdelBio2+0.5*dzdelBio2)/dzdelBio2;
        end
    end

    %Biological Parameters
    ws2 = 0.5/86400; %[m/s] Phytoplankton sinking rate. Negative built in code.
    a2 = 0.1*86400; %[m^2*s/mol quanta] Photosynthetic efficiency at low irradiance
    I02 = 30/86400; %[mol quanta m^2/s] Mean daily surface irradiance
    Pmax2 = 100/86400; %[mgC(mg cla)^-1 s^-1] Maximum phytoplankton carbon assimilation rate
    r2 = 0.05; %[% of Pmax] Respiration rate 
    ZP2 = 0.1/86400; %[s^-1] Zooplankton grazing rate
    BG2 = 0.1/86400; %Fix for now. %[m^3 m^-2 s^-1] Benthic grazing rate 0.1
    theta2 = 50; %[mgC(mg chl a)^-1] Ratio phytoplankton cellular carbon to chlorophyll
    kt2 = 2; %[m-1] Abiotic light attenuation coefficient. Range 1-7. Used to be 3, now 2
    if want_self_shade == 1
        kb2_coeff = 0.016; %Coefficient for self-shading
    else
        kb2_coeff = 0;
    end
    if irr_type == 2
        kb2 = kb2_coeff*cumsum(Bio2); %[m-1] Self-shading by phyto-plankton biomass (correct way)
    else
        kb2 = kb2_coeff*Bio2; %[m-1] Self-shading by phyto-plankton biomass
    end
    
    if irradiance == 1 %If want to perform irradiance
        if irr_type == 3
            %I2 = I02.*exp((-kt2-kb2).*abs(z2));
            I2(1) = I02*exp(-(kt2+kb2(1))*dz/2);
            for i=2:N2
                I2(i)=I2(i-1)*exp(-(kt2+kb2(i))*dz);
            end
        else
            I2 = I02.*exp((-kt2).*abs(z2)-kb2);
        end
        I2_flipped = flip(I2); %Flip so high irradiance at top
       % f_I2 = tanh(a2*I2);
        f_I2 = tanh(a2*I2_flipped);
    elseif irradiance == 2 %If want to perform Roy irradiance
        %Irradiance - parameters from Roy
        Vol2 = 1.6; %Average body volume phytoplankton
        Ec2 = 0.2; %Extinction coefficient water
        Ext2 = 0.12*Vol2^-0.33; %Self shading of phytoplankton due to body size
        Iopt2 = 680; %Optimum surface solar radiation for photosynthesis
        Ir2 = 800; %Surface solar irradiance.

        %Roy irradiance effect
        I2 = Ir2.*exp((-Ec2-Ext2)*abs(z2)); 
        
        f_I2 = I2./Iopt2.*exp(1-I2./Iopt2); 
    else
        f_I2 = ones(1,N2);
    end
    
    munet2 = Pmax2*(f_I2-r2)/theta2;
    
    % *******************************************************************
    %  Set initial conditions for u, q^2, and other turbulence quantities
    % *******************************************************************
    for i=1:N2
       U2(i) = 0.0;
       V2(i) = 0.0;
       Q22(i)=SMALL;  %"seed" the turbulent field with small values, then let it evolve
       Q2L2(i)=SMALL;
       Q_2(i)=sqrt(Q22(i));
       L2(i)=-kappa*H2*(z2(i)/H2)*(1-(z2(i)/H2)); %Q2L(n,1)/Q2(n,1) = 1 at initialization;
       Gh2= -(N_BV2(i)*L2(i)/(Q_2(i)+SMALL))^2; 
       Gh2=min(Gh2,0.0233);
       Gh2=max(Gh2,-0.28);
       num2=B12^(-1/3)-A12*A22*Gh2*((B22-3*A22)*(1-6*A12/B12)-3*C12*(B22+6*A12));
       dem2=(1-3*A22*Gh2*(B22+6*A12))*(1-9*A12*A22*Gh2);
       Sm2(i)=num2/dem2;
       Sh2(i)=A22*(1-6*A12/B12)/(1-3*A22*Gh2*(B22+6*A12)); 
       Kq2(i)=Sq2*Q_2(i)*L2(i)+nu;  %turbulent diffusivity for Q2
       nu_t2(i)=Sm2(i)*Q_2(i)*L2(i)+nu; %turbulent viscosity
       Kz2(i)=Sh2(i)*Q_2(i)*L2(i)+ nu; %turbulent scalar diffusivity
    end

    % *******************************************************************
    %  Save initial conditions as first columns in saved matrix
    % *******************************************************************
    Um2(:,1) = U2; %Velocity
    Cm2(:,1) = C2; %Salinity
    Biom2(:,1) = Bio2; %Biomass
    Q2m2(:,1) = Q22; %TKE
    Q2Lm2(:,1) = Q2L2; %TKE*Length scale
    rhom2(:,1) = rho2; %Density
    Lm2(:,1) = L2; %Length scale
    nu_tm2(:,1) = nu_t2; %Eddy viscosity
    Kzm2(:,1) = Kz2; %Eddy diffusivity
    Kqm2(:,1) = Kq2; %Effective diffusion coefficient TKE
    N_BVm2(:,1) = N_BV2; %Brunt-Vaisala frequency
end

% *******************************************************************
%%  Set initial conditions for Ky and KyS- outside flag so still runs
% *******************************************************************
dy = 1000; %[m] Width of channel
dy2 = 4000; %[m] Width of shoal
Ly = 0.5*(dy+dy2); %[m] Distance between center WC1 and WC2
for i=1:N2
    if imposeKy == 0
        Ky(i) = Ky_range(Kytest); %[m^2/s] Effective rate of Bio diffusive exchange 
        KyS(i) = Ky_range(Kytest); %Effective rate of salinity exchange
    else %Use average Ky value to compare
        Ky(i) = Kyave;
        KyS(i) = Kyave;
    end
end

Kymin = Kymin_range(Kytest); %Take minimum value for cycling through.

%Save information about Ky without needing whole vector
Kytop = Ky(N2);
Kybottom = Ky(1);
Kytype = 'constant';

KyStop = KyS(N2);
KySbottom = KyS(1);
KyStype = 'constant';

ratio1 = 1/(Ly*dy); %Ratio of concentration exchange from WC2 to WC1. [1/m^2]
ratio2 = 1/(Ly*dy2); %Ratio of concentration exchange from WC1 to WC2

Kym(:,1) = Ky; 
KySm(:,1) = KyS;
Kyp = Ky;
KySp = KyS;

Umax = 0.5; %[m/s] Maximum velocity value to permit in condition 3 below.
num_range = [0,0.8750,1.75,3.5,7,14,28,56,12.4/24,1]; %First value to show that didn't use this variable this run [0,7,14,28,56];  [0,0.8750,1.75,3.5,7,14,28,56]
ebb_range = [-1,0,1]; %-1 means didn't use this variable this run, 0 means ebb, 1 means flood
numdays = num_range(daystest); %Number days for Ky cycle condition 2. Works by turning off numdays, then on numdays, etc.
ebbfl = ebb_range(eftest);

%Conditions for Ky time variable. Make in line function with @(t,U) even if
%no time or U dependence so will run in wc_setup. Make a logical.
if condnum == 1 
    if ebbfl == 0
        condition = @(t,U) sign(mean(U)) < 0; %Ebb is when < 0
    elseif ebbfl == 1
        condition = @(t,U) sign(mean(U)) > 0; %Flood is when > 0
    end
elseif condnum == 2 
    if ebbfl == 0 || ebbfl == -1 %In sync with spring neap (natural). Default case
        condition = @(t,U) sign(cos(2*pi*t/(3600.*24.*numdays.*2))) > 0; 
    elseif ebbfl == 1 %Out of sync with spring neap (gates)
        condition = @(t,U) sign(cos(2*pi*t/(3600.*24.*numdays.*2))) < 0; 
    end
elseif condnum == 3
    condition = @(t,U) abs(mean(U)) < Umax; %Slack tide
elseif condnum == 4 
    if ebbfl == 0
        condition = @(t,U) sign(mean(U)) < 0 && abs(mean(U)) < Umax; %Ebb and slack
    elseif ebbfl == 1 
        condition = @(t,U) sign(mean(U)) > 0 && abs(mean(U)) < Umax; %Flood and slack
    end
end

% *******************************************************************
%%  Set up calculate difference
% *******************************************************************
if first_water_column==1 && second_water_column==1
    diffUm(:,1) = U2 - U(N-N2+1:N); %Top section of WC1 is shoal
    diffCm(:,1) = C2 - C(N-N2+1:N);
    diffBiom(:,1) = Bio2 - Bio(N-N2+1:N);
    diffQ2m(:,1) = Q22 - Q2(N-N2+1:N);
    diffQ2Lm(:,1) = Q2L2 - Q2L(N-N2+1:N);
    diffrhom(:,1) = rho2 - rho(N-N2+1:N);
    diffLm(:,1) = L2 - L(N-N2+1:N);
    diffnu_tm(:,1) = nu_t2 - nu_t(N-N2+1:N);
    diffKzm(:,1) = Kz2 - Kz(N-N2+1:N);
    diffKqm(:,1) = Kq2 - Kq(N-N2+1:N);
    diffN_BVm(:,1) = N_BV2 - N_BV(N-N2+1:N);
    
    %Exchange
    Y = [diffBiom(:,1)';diffCm(:,1)'];
    for j=1:2  %Bio,C
        if j < 2
            exchange(j,:) = Ky.*Y(j,:); %Bio exchange
        else
            exchange(j,:) = KyS.*Y(j,:); %Salinity exchange
        end
    end
end

%CHANGE TESTS: FIXME
coeff = 0.01; %0.025
coeff2 = 0.19;
Cx = 0.00025; %0.00025
Cx2 = 0.00025;
Px0 = 1.25*Px0;
Px02 = 1.25*Px02;

% ws = 0.5/86400; %[m/s]
% ws2 = 0.5/86400; %[m/s]
Px0 = 1.25*Px0;
Px02 = 1.25*Px02; 

% kt2 = 2; %[m-1] Abiotic light attenuation coefficient. Range 1-7 Low 2, Mid 2.5

% *******************************************************************
%%  Save initial conditions to .MAT file
% *******************************************************************
num_sets = 0; %Which set of tests this is on
run_num = 0; %Which run of this set this is on

if saveme == 1  && loadme == 0 %Whether or not to save this run. Don't save if loading file
    if isfile([savefile,'\',filename,'.mat']) %See if this file exists
        if allvars == 1
            save('lastrun')      %Save initial variables this run
        else
            for i = 1:numel(savewhat)
                if i == 1
                    save('lastrun', savewhat{i}); %Save certain variables
                else
                    save('lastrun', savewhat{i},'-append'); %Save certain variables
                end
            end 
        end
        
        thisrun = load('lastrun'); %Save as matrix
        mmm = matfile([savefile,'\',filename,'.mat']); %Evaluate details of file without loading: save time

        if isprop(mmm,'vars')==1  %By third run, vars exist. Can just load 
            load([savefile,'\',filename,'.mat']); 
        else
            vars = load([savefile,'\',filename,'.mat']); %Create struct of vars
        end
        
        vars = [vars,thisrun];
        [num_sets,run_num] = size(vars);

        save([savefile,'\',filename,'.mat'],'vars'); %Save new variable set

        if export == 1 %Export struct to excel
            writetable(struct2table(vars),[savefile,'\',filename,'.xlsx']); 
        end
        delete('lastrun.mat')
    else %First run save new file and directory
        if isfolder(filename) == 0 %Make directory to save plot outputs
            mkdir(filename) 
        end
        if allvars == 1
            save([savefile,'\',filename,'.mat'])     %Save initial variables this run
        else
            for i = 1:numel(savewhat)
                if i == 1
                    save([savefile,'\',filename,'.mat'], savewhat{i}); %Save certain variables
                else
                    save([savefile,'\',filename,'.mat'], savewhat{i},'-append'); %Save certain variables
                end
            end 
        end
        num_sets = 1; %Initialize 
        run_num = 1;
    end
    disp(['File saved. File: ', filename, '. Run_num: ' num2str(run_num)])
end

% *******************************************************************
%%  Load initial conditions from a .MAT file (replace what just made)
% *******************************************************************
if loadme == 1 && saveme == 0 %Whether or not to take initial conditions from loaded file. Cannot be loaded if saved.
    mmm = matfile([savefile,'\',filename,'.mat']); %Evaluate details of file without loading: save time
    
    if isprop(mmm,'vars')==1  %By third run, vars exist. Can just load 
        load([savefile,'\',filename,'.mat']); 
    else
        vars = load([savefile,'\',filename,'.mat']); %Create struct of vars
    end
    
    names = fieldnames(vars);
    run_num = loadnum;
    for i = 1:numel(names)
        assignin('caller', names{i}, vars(loadnum).(names{i}));
    end
    
    %Recreate Ky from loaded file 
    if strcmp(Kytype,'constant') == 1
        if imposeKy == 1
            Ky = ones(1,N2)*Kyave; %Impose average Ky value from previous run
            KyS = ones(1,N2)*Kyave;
            Kyp = Ky;
            KySp = KyS;
        else
            Ky = ones(1,N2)*Kytop;
            KyS = ones(1,N2)*KyStop;
            Kyp = Ky;
            KySp = KyS;
        end
    else
        disp(['Ky type not defined. Treated as ',num2str(Ky_range(Kytest))])
    end
    
    %Recreate Ky time variable
    if condnum == 1 
        if ebbfl == 0
            condition = @(t,U) sign(mean(U)) < 0; %Ebb is when < 0
        elseif ebbfl == 1
            condition = @(t,U) sign(mean(U)) > 0; %Flood is when > 0
        end
    elseif condnum == 2
        if ebbfl == 0 || ebbfl == -1 %In sync with spring neap (natural). Default case
            condition = @(t,U) sign(cos(2*pi*t/(3600.*24.*numdays.*2))) > 0; 
        elseif ebbfl == 1 %Out of sync with spring neap (gates)
            condition = @(t,U) sign(cos(2*pi*t/(3600.*24.*numdays.*2))) < 0; 
        end
    elseif condnum == 3
        condition = @(t,U) abs(mean(U)) < Umax; %Slack tide
    elseif condnum == 4 
        if ebbfl == 0
            condition = @(t,U) sign(mean(U)) < 0 && abs(mean(U)) < Umax; %Ebb and slack
        elseif ebbfl == 1 
            condition = @(t,U) sign(mean(U)) > 0 && abs(mean(U)) < Umax; %Flood and slack
        end
    end
    
    if mod(loadnum,2)==0 %Recreating averaged Ky case
        imposeKy = 1;
    end
    
    %FIXME: need to recreate other changed parameters such as P0 -> initial 
    %P vector. Load beginning or dialog box? 
    
    disp(['File loaded. File: ', filename, '. Run_num: ' num2str(loadnum)])
end

%CHANGE TESTS AFTER LOADFILE
% ws = 0.5/86400; %[m/s]
% ws2 = 0.5/86400; %[m/s]
% Px0 = 1.25*Px0;
% Px02 = 1.25*Px02; 
% 
% kt2 = 2.5; %[m-1] Abiotic light attenuation coefficient. Range 1-7

%Test sinking flux: 
% Px0 = 0;
% Px02 = 0;
% BG = 0;
% BG2 = 0; 
% Pmax = 0;
% Pmax2 = 0;
% munet = Pmax*(f_I2-r2)/theta2;
% munet2 = Pmax2*(f_I2-r2)/theta2;
% ZP = 0;
% ZP2 = 0;
% ws = 0;
% ws2 = 0;

%condition = @(t,U) sign(cos(2*pi*t/(3600.*24.*numdays.*2))) < 0; 