% ********************************************************************
% Master code for 200B water column (vertical) cases
%    Original code from Lisa Lucas, modified by Tina Chow 
%       Spring 2018: Mark Stacey and Michaella Chung
%       Summer & Fall 2019 modifications: Lily Engel
%       Now needs to be initiated from file wc_runme_developmodelv3.m
% ********************************************************************
close all;

% ********************************************************************
%% Define model set up - grid and timestep
% ********************************************************************
N=80;%number of grid points 80
H=15; %depth WC1-channel (meters)
H2=2; %depth WC2-shoal (meters)
dz=H/N; %grid spacing - may need to adjust to reduce oscillations. Use same for both WCs so line up.
N2=round(H2/dz); %number of grid points WC2 to keep same dz 
dt=60; %(seconds) size of time step 
if daystest == 8 %224 for 56 days extended test
    M=224*24*60; %number of time steps 
else
    M=112*24*60; %number of time steps 3600, 1440, 2400, 31205 (*2), 20160, 72000, (84: 120960, 112, 224)
    %M=224*24*60; %number of time steps 3600, 1440, 2400, 31205 (*2), 20160, 72000, (84: 120960, 112, 224)
end
beta=dt/dz^2;
z = zeros(1,N); %Pre-allocate for speed. WC1
z2 = zeros(1,N2); %Pre-allocate for speed. WC2

for i=1:N % Initialize grid   
   z(i)=-H+dz*(i-1/2); %bottom at z=-H, free surface at 0
   if i<=N2
       z2(i)=-H2+dz*(i-1/2); %bottom at z2=-H2, free surface at 0
   end
end
isave=1; %increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
savecount=0;

first_water_column=1; %1 means run, anything else means not
second_water_column=1;

% *******************************************************************
%%  Call code to set up initial conditions and model parameters
% *******************************************************************
wc_preallocate %Pre-allocate large matrices for speed
wc_setup

if loadmePhy ~= 1 || loadme ~= 0 %Only run time loop if not recreating something already from a file
    % *******************************************************************
    %  Start of time loop
    % *******************************************************************
    for m=2:M
        t(m)=dt*(m-1); %define time
        wc_advance %uses BGO/Mellor-Yamada 2-equation closure (level 2.5)

    % *******************************************************************
    %  Saving Profiles - isave defines decimation
    % ******************************************************************
        if mod(m,isave) == 0
            savecount = savecount+1;
            times(savecount) = t(m); %#ok<*SAGROW> 
            
            %FIRST WATER COLUMN 
            if first_water_column == 1
                Um(:,savecount) = U; %Velocity
                Cm(:,savecount) = C; %Salinity
                Biom(:,savecount) = Bio; %Biomass
                Q2m(:,savecount) = Q2; %TKE
                Q2Lm(:,savecount) = Q2L; %TKE*Length scale
                rhom(:,savecount) = rho; %Density
                Lm(:,savecount) = L; %Length scale
                nu_tm(:,savecount) = nu_t; %Eddy viscosity
                Kzm(:,savecount) = Kz; %Eddy diffusivity
                Kqm(:,savecount) = Kq; %Effective Diffusion coefficient
                N_BVm(:,savecount) = N_BV; %Brunt-Vaisala frequency
                Cplot(savecount) = Cm(N,savecount)-Cm(1,savecount);
                Pxm(:,savecount) = Px; %Pressure-forcing term
                munetm(:,savecount) = munet;
            end

            %SECOND WATER COLUMN
            if second_water_column == 1
                Um2(:,savecount) = U2; 
                Cm2(:,savecount) = C2; 
                Biom2(:,savecount) = Bio2;
                Q2m2(:,savecount) = Q22;
                Q2Lm2(:,savecount) = Q2L2;
                rhom2(:,savecount) = rho2;
                Lm2(:,savecount) = L2; 
                nu_tm2(:,savecount) = nu_t2; 
                Kzm2(:,savecount) = Kz2; 
                Kqm2(:,savecount) = Kq2; 
                N_BVm2(:,savecount) = N_BV2; 
                Cplot2(savecount) = Cm2(N2,savecount)-Cm2(1,savecount);
                Pxm2(:,savecount) = Px2;
                munetm2(:,savecount) = munet2;
            end

            %CALCULATE DIFFERENCE
            if first_water_column == 1 && second_water_column == 1
                Kym(:,savecount) = Ky;
                KySm(:,savecount) = KyS;
                diffUm(:,savecount) = diffU;
                diffCm(:,savecount) = diffC;
                diffBiom(:,savecount) = diffBio;
                diffQ2m(:,savecount) = diffQ2;
                diffQ2Lm(:,savecount) = diffQ2L;
                diffrhom(:,savecount) = diffrho;
                diffLm(:,savecount) = diffL;
                diffnu_tm(:,savecount) = diffnu_t;
                diffKzm(:,savecount) = diffKz;
                diffKqm(:,savecount) = diffKq;
                diffN_BVm(:,savecount) = diffN_BV;
            end
            
            %Make Movie
            if movie == 1 && mod(m,1000) == 0
                fig = figure(41);
                timestr = sprintf('%5.3f',t(m)/86400); %Convert the current time to a string for printing in the title
                subtightplot(4,4,[1 13])
                BB = Bio'.*ones(N);
                contourf(y1,z,BB,'Linestyle','none')
                if wantsamecolor == 1 
                    caxis manual
                    caxis([bottom top]);
                end
                colormap(colorscheme)
                ylabel('Depth (m)')
                xlabel('Width (m)')
                if Ky(1) == Kymin
                    title('Channel Exchange OFF')
                else
                    title('Channel Exchange ON')
                end
                subtightplot(4,4,[2 16])
                BB2 = Bio2'.*ones(N2);
                contourf(y2,z2,BB2,'Linestyle','none')
                colormap(colorscheme)
                ylim([z(1) z2(N2)])
                set(gca,'YTick',[])
                xlabel('Width (m)')
                c = colorbar('south');
                c.Label.String = 'Biomass Conc (mg chl a/m^3)';
                title(['Shoal Case Study t = ',timestr,'days'])
                if wantsamecolor == 1 
                    caxis manual
                    caxis([bottom top]);
                end
                set(fig,'DoubleBuffer','on'); %this speeds up the display of the image
                set(gca,'nextplot','replace','Visible','on') %Improves display and capture of the image in the figure window.
                %sgtitle([filename, ' Run ', num2str(run_num),' t = ',timestr,'d']); %Display this title
                frame = getframe(gcf); %Capture the frame
                writeVideo(aviobj1,frame); %Add the frame to your avi movie file
                clf
            end
        end
    end

    % *******************************************************************
    %  End of time loop
    % *******************************************************************
end

if movie==1
    close(aviobj1);
end

% *******************************************************************
%%  PLOTTING AND CALCULATING FOLLOWS HERE
%   columns of variable matrices (Um, Cm, etc) vs. z array
% *******************************************************************
%Calculations
if first_water_column == 1
    if si_type == 1
        Si = g*abs(alpha)*Cx*H^2/(C_D*(1/(2*pi)*Px0*3600*T_Px)^2); %Simpson Number new way - preferred
    else
        Si_opt2 = g*abs(alpha)*Cx*H^2/(C_D*max(max(abs(Um)))^2); %Simpson number Monismith et. al (PS 4 w/C_D - includes baroclinic term already in Um)
    end

    Un = 2*pi*H/(T_Px*3600*sqrt(C_D)*(1/(2*pi)*Px0*3600*T_Px)); %Unsteadiness Number

    if calc_rous == 1
        Ro = ws/(kappa*sqrt(C_D)*(1/(2*pi)*Px0*3600*T_Px)); %Rouse Number
    end
end

if second_water_column == 1
    if si_type == 1
        Si2 = g*abs(alpha2)*Cx2*H2^2/(C_D*(1/(2*pi)*Px02*3600*T_Px2)^2); %Simpson Number new way - preferred
    else
        Si_opt22 = g*abs(alpha2)*Cx2*H2^2/(C_D*max(max(abs(Um2)))^2); %Simpson number Monismith et. al (PS 4 w/C_D - includes baroclinic term already in Um)
    end

    Un2 = 2*pi*H2/(T_Px2*3600*sqrt(C_D)*(1/(2*pi)*Px02*3600*T_Px2)); %Unsteadiness Number

    if calc_rous == 1
        Ro2 = ws2/(kappa*sqrt(C_D)*(1/(2*pi)*Px02*3600*T_Px2)); %Rouse Number
    end
end

if first_water_column == 1 && second_water_column ==1 
    %Depth averaged Bio
    depaBio = 1/H*sum(dz.*Biom);
    depaBio2 = 1/H2*sum(dz.*Biom2);
    if print == 1  %Print all the following 
        fprintf(['Alpha = ', num2str(alpha),',',num2str(alpha2),'\n'])
        fprintf(['\nT_Px = ', num2str(T_Px),',',num2str(T_Px2),'\n'])
        fprintf(['\nPx = ', num2str(Px0),',',num2str(Px02),'\n'])
        fprintf(['\nCx = ', num2str(Cx),',',num2str(Cx2),'\n'])
        if si_type==1 
            fprintf(['\nSi = ', num2str(Si),',',num2str(Si2),'\n'])
        else
            fprintf(['\nSi_opt2 = ', num2str(Si_opt2),',',num2str(Si_opt22),'\n'])
        end
        fprintf(['\nUn = ', num2str(Un),',',num2str(Un2),'\n'])
        if calc_rous == 1
            fprintf(['\nRo = ', num2str(Ro),',',num2str(Ro2),'\n'])
        end
        fprintf(['\nKy = ', num2str(Ky_range(Kytest)),'\n'])
    end
    bottom = min(min(min(Biom)),min(min(Biom2))); %Min and max of bio colorbars
    top  = max(max(max(Biom)),max(max(Biom2)));
end

%Plots
first = times(1)/3600; 
last = times(length(times))/3600; %If want to plot whole time, do first:last

tstart = find(times/3600==first); %First hour want to plot. 100, 42*24, first, 56*24
tend = find(times/3600==last); %Last hour want to plot. last, 70*24

if daystest == 8 %224 for 56 days extended test
    teqstart = find(times/3600==112*24); %First hour for equilibrium want to plot 42, 56 (112 for 56 days test)
else
    teqstart = find(times/3600==56*24); %First hour for equilibrium want to plot 42, 56 
 %   teqstart = find(times/3600==112*24); %First hour for equilibrium want to plot 42, 56 
end

if saveme == 0 && loadme == 0 && loadmePhy == 0 %If hardcoding to make plots
    run_num = 0; %Whichever run recreating without loading or saving
end

whichplot = 1; %Run first group of plots
plots %Run plotting code

% *******************************************************************
%  End of Main Program
% *******************************************************************