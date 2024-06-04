% *******************************************************************
%  PLOTTING FOLLOWS HERE
%   Three different sections for different parts of the code. All operated
%   by if statements FIXME - switch to cases? Or make a function?
% *******************************************************************

%% First set plots
if whichplot == 1
    %FIRST WATER COLUMN-------------------------------------------------------
    if first_water_column == 1
        if plotme(1) == 1 %WC1 Salinity
            figure(1)
            contourf(times(tstart:tend)/3600,z,Cm(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title({[filename, ' Run ', num2str(run_num)];['WC 1 Salinity Cx = ',num2str(Cx),', coeff= ',num2str(coeff)]})
            c = colorbar;
            c.Label.String = 'Salinity (psu)';
        end

        if plotme(2) == 1 %WC1 Velocity
            figure(2)
            contourf(times(tstart:tend)/3600,z,Um(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title({[filename, ' Run ', num2str(run_num)];'WC 1 Velocity'})
            c = colorbar;
            c.Label.String = 'Velocity (m/s)';
        end

        if plotme(3) == 1 %WC1 TKE
            figure(3)
            contourf(times(tstart:tend)/3600,z,Q2m(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title({[filename, ' Run ', num2str(run_num)];'WC 1 Turbulent Kinetic Energy'})
            c = colorbar;
            c.Label.String = 'TKE (m^2/s^2)'; 
        end

        if plotme(4) == 1 %WC1 Bio surf
            figure(4)
            surf(times/3600,z,Biom)
            shading interp
            colormap(colorscheme)
            c = colorbar;
            c.Label.String = 'Biomass Conc (mg chl a/m^3)'; 
            xlabel('Time (hrs)')
            ylabel('Depth (m)')
            zlabel('Biomass Conc (mg chl a/m^3)')
            title({[filename, ' Run ', num2str(run_num)];'WC 1 Biomass'})
        end

        if plotme(5) == 1 %WC1 Bio sub
            figure(5)
            contourf(times(tstart:tend)/3600,z,Biom(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title({[filename, ' Run ', num2str(run_num)],['WC 1 Biomass, Biox = ',num2str(Biox)]})
            colorbar
            c = colorbar;
            c.Label.String = 'Biomass Conc (mg chl a/m^3)'; 
        end

        if plotme(19) == 1 %WC1 Bottom Boundary
            for j=tstart:tend
                figure(19)
           %     j=tstart; %Time to plot
                plot(Biom(:,j),z)
                ylabel('Depth (m)')
                xlabel('Biomass (mg chl a/m^3)')
                title({[filename, ' Run ', num2str(run_num)],['WC 1 Bottom Boundary Time step = ',num2str(j)]})
            end
        end

        if plotme(21) == 1 %WC1 Eddy Diffusivity
            figure(21)
            contourf(times(tstart:tend)/3600,z,log(Kzm(:,tstart:tend)),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title({[filename, ' Run ', num2str(run_num)];['Log WC 1 Kz Cx = ',num2str(Cx),' coeff= ',num2str(coeff)]})
            c = colorbar;
            c.Label.String = 'Log(Kz) (m^2/s)'; 
        end
    end

    %SECOND WATER COLUMN-------------------------------------------------------
    if second_water_column == 1
        if plotme(6) == 1 %WC2 Salinity
            figure(6)
            contourf(times(tstart:tend)/3600,z2,Cm2(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title({[filename, ' Run ', num2str(run_num)];['WC 2 Salinity Cx = ',num2str(Cx2),', coeff= ',num2str(coeff2)]})
            c = colorbar;
            c.Label.String = 'Salinity (psu)'; 
        end

        if plotme(7) == 1 %WC2 Velocity
            figure(7)
            contourf(times(tstart:tend)/3600,z2,Um2(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title({[filename, ' Run ', num2str(run_num)];'WC 2 Velocity'})
            c = colorbar;
            c.Label.String = 'Velocity (m/s)'; 
        end

        if plotme(8) == 1 %WC2 TKE
            figure(8)
            contourf(times(tstart:tend)/3600,z2,Q2m2(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title({[filename, ' Run ', num2str(run_num)];'WC 2 Turbulent Kinetic Energy'})
            c = colorbar;
            c.Label.String = 'TKE (m^2/s^2)'; 
        end

        if plotme(9) == 1 %WC2 Bio surf
            figure(9)
            surf(times/3600,z2,Biom2)
            shading interp
            colormap(colorscheme)
            c = colorbar;
            c.Label.String = 'Biomass Conc (mg chl a/m^3)'; 
            xlabel('Time (hrs)')
            ylabel('Depth (m)')
            zlabel('Biomass Conc (mg chl a/m^3)')
            title({[filename, ' Run ', num2str(run_num)];'WC 2 Biomass'})
        end

        if plotme(10) == 1 %WC2 Bio sub
            figure(10)
            contourf(times(tstart:tend)/3600,z2,Biom2(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title({[filename, ' Run ', num2str(run_num)];['WC 2 Biomass, Biox = ',num2str(Biox2)]})
            c = colorbar;
            c.Label.String = 'Biomass Conc (mg chl a/m^3)'; 
        end

        if plotme(20) == 1 %WC2 Bottom Boundary
            for j=tstart:tend
                figure(20)
               % j=tstart; %Time to plot
                plot(Biom2(:,j),z2)
                ylabel('Depth (m)')
                xlabel('Biomass (mg chl a/m^3)')
                title({[filename, ' Run ', num2str(run_num)],['WC 2 Bottom Boundary Time step = ',num2str(j)]})
            end
        end

        if plotme(22) == 1 %WC2 Eddy Diffusivity
            figure(22)
            contourf(times(tstart:tend)/3600,z2,log(Kzm2(:,tstart:tend)),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title({[filename, ' Run ', num2str(run_num)];['Log WC 2 Kz Cx = ',num2str(Cx2),' coeff= ',num2str(coeff2)]})
            c = colorbar;
            c.Label.String = 'Log(Kz) (m^2/s)'; 
        end
    end

    %BOTH WATER COLUMNS--------------------------------------------------------
    if first_water_column == 1 && second_water_column == 1
        if plotme(11) == 1 %Both WCs Salinity
            figure(11)
            sgtitle([filename, ' Run ', num2str(run_num)])
            subplot(1,2,1)
            contourf(times(tstart:tend)/3600/24,z,Cm(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (days)')
            title(['WC 1 Salinity Cx = ',num2str(Cx),', coeff= ',num2str(coeff)])
            c = colorbar;
            c.Label.String = 'Salinity (psu)'; 
            hold on
            yyaxis right
            plot(times(tstart:tend)/3600/24,depaBio(tstart:tend),'k-','LineWidth',3)
            ylabel('Depth- Averaged Biomass') 
            hold off
            subplot(1,2,2)
            contourf(times(tstart:tend)/3600/24,z2,Cm2(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
          %  ylim([z(1) z2(N2)])
            xlabel('Time (days)')
            title(['WC 2 Salinity Cx = ',num2str(Cx2),', coeff= ',num2str(coeff2)])
            c = colorbar;
            c.Label.String = 'Salinity (psu)'; 
            hold on
            yyaxis right
            plot(times(tstart:tend)/3600/24,depaBio2(tstart:tend),'k-','LineWidth',3)
            ylabel('Depth- Averaged Biomass') 
            hold off
            
            figure()
            sgtitle([filename, ' Run ', num2str(run_num)])
            subplot(1,2,1)
            contourf(times(tstart:tend)/3600/24,z(2:N),diff(Cm(:,tstart:tend))/dz,'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (days)')
            title(['WC 1 dSdz Cx = ',num2str(Cx),', coeff= ',num2str(coeff)])
            c = colorbar;
            c.Label.String = 'Salinity (psu)'; 
            hold on
            yyaxis right
            plot(times(tstart:tend)/3600/24,depaBio(tstart:tend),'k-','LineWidth',3)
            ylabel('Depth- Averaged Biomass') 
            hold off
            subplot(1,2,2)
            contourf(times(tstart:tend)/3600/24,z2(2:N2),diff(Cm2(:,tstart:tend))/dz,'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
          %  ylim([z(1) z2(N2)])
            xlabel('Time (days)')
            title(['WC 2 dSdz Cx = ',num2str(Cx2),', coeff= ',num2str(coeff2)])
            c = colorbar;
            c.Label.String = 'Salinity (psu)'; 
            hold on
            yyaxis right
            plot(times(tstart:tend)/3600/24,depaBio2(tstart:tend),'k-','LineWidth',3)
            ylabel('Depth- Averaged Biomass') 
            hold off
        end

        if plotme(12) == 1 %Both WCs Velocity
            figure(12)
            sgtitle([filename, ' Run ', num2str(run_num)])
            subplot(1,2,1)
            contourf(times(tstart:tend)/3600,z,Um(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title('WC 1 Velocity')
            c = colorbar;
            c.Label.String = 'Velocity (m/s)'; 
            subplot(1,2,2)
            contourf(times(tstart:tend)/3600,z2,Um2(:,tstart:tend),'Linestyle','none');
            colormap(colorscheme)
            ylabel('Depth (m)')
            ylim([z(1) z2(N2)])
            xlabel('Time (hrs)')
            title('WC 2 Velocity')
            c = colorbar;
            c.Label.String = 'Velocity (m/s)'; 
        end

        if plotme(13) == 1 %Both WCs TKE
            figure(13)
            sgtitle([filename, ' Run ', num2str(run_num)])
            subplot(1,2,1)
            contourf(times(tstart:tend)/3600,z,Q2m(:,tstart:tend));
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (hrs)')
            title('WC 1 Turbulent Kinetic Energy')
            c = colorbar;
            c.Label.String = 'TKE (m^2/s^2)'; 
            subplot(1,2,2)
            contourf(times(tstart:tend)/3600,z2,Q2m2(:,tstart:tend));
            colormap(colorscheme)
            ylabel('Depth (m)')
            ylim([z(1) z2(N2)])
            xlabel('Time (hrs)')
            title('WC 2 Turbulent Kinetic Energy')
            c = colorbar;
            c.Label.String = 'TKE (m^2/s^2)'; 
        end
 
        if plotme(14) == 1 %Both WCs Bio surfplots
            figure(14)
            sgtitle([filename, ' Run ', num2str(run_num)])
            subplot(1,2,1)
            surf(times/3600,z,Biom)
            shading interp
            colormap(colorscheme)
            if wantsamecolor == 1 
                caxis manual
                caxis([bottom top]);
            end
            c = colorbar;
            c.Label.String = 'Biomass Conc (mg chl a/m^3)'; 
            xlabel('Time (hrs)')
            ylabel('Depth (m)')
            zlabel('Biomass Conc (mg chl a/m^3)')
            title('WC 1 Biomass')
            subplot(1,2,2)
            surf(times/3600,z2,Biom2)
            shading interp
            colormap(colorscheme)
            if wantsamecolor == 1 
                caxis manual
                caxis([bottom top]);
            end
            c = colorbar;
            c.Label.String = 'Biomass Conc (mg chl a/m^3)';
            xlabel('Time (hrs)')
            ylabel('Depth (m)')
            ylim([z(1) z2(N2)])
            zlabel('Biomass Conc (mg chl a/m^3)')
            title('WC 2 Biomass')
        end

        if plotme(15) == 1 %Both WCs Bio subplots
            figure(15)
            sgtitle([filename, ' Run ', num2str(run_num)])
            subplot(1,2,1)
            contourf(times(tstart:tend)/3600/24,z,Biom(:,tstart:tend),'Linestyle','none')
            if wantsamecolor == 1 
                caxis manual
                caxis([bottom top]);
            end
            colormap(colorscheme)
            ylabel('Depth (m)')
            xlabel('Time (Days)')
            title('WC 1 Biomass')
            c = colorbar;
            c.Label.String = 'Biomass Conc (mg chl a/m^3)';
            subplot(1,2,2)
            contourf(times(tstart:tend)/3600/24,z2,Biom2(:,tstart:tend),'Linestyle','none')
            colormap(colorscheme)
            ylabel('Depth (m)')
            ylim([z(1) z2(N2)])
            xlabel('Time (Days)')
            title('WC 2 Biomass')
            if wantsamecolor == 1 
                caxis manual
                caxis([bottom top]);
            end
            c = colorbar;
            c.Label.String = 'Biomass Conc (mg chl a/m^3)';
        end

        if plotme(16) == 1 %Both WCs Middle water column
            figure(16)
            dep = -4;
            j = find(z >= dep, 1); %Takes the first of what indices it finds
            if abs(dep) < H2
                plot(times(tstart:tend)/3600,Biom(j,tstart:tend),times(tstart:tend)/3600,Biom2(j-N+N2,tstart:tend))
                legend('B1','B2')
                title({[filename, ' Run ',num2str(run_num)]; ['Fixed Depth z = ',num2str(dep)]})
                ylabel('Biomass (mg chl a/m^3)')
                xlabel('Time (hrs)')
            else
                disp('Depth not in shoal.')
            end
        end

        if plotme(17) == 1 %Both WCs Depth Averaged Bio
            figure(17)
            subplot(1,2,1)
            plot(times(tstart:tend)/3600/24,depaBio(tstart:tend),'DisplayName','Tvar Ky')
            ylabel('Biomass (mg chl a/m^3)')
            xlabel('Time (Days)')
            title('Channel')
            subplot(1,2,2)
            plot(times(tstart:tend)/3600/24,depaBio2(tstart:tend),'r-','DisplayName','Tvar Ky')
            sgtitle({[filename, ' Run ',num2str(run_num)]; ['Depth Averaged Bio, Max Ky = ',num2str(Kytop), ' m^2/s']})
            ylabel('Biomass (mg chl a/m^3)')
            xlabel('Time (Days)')
            title('Shoal')
        end

        if plotme(18) == 1 %Separate fluxes of Bio expression. 
            f_I = 1/H*sum(dz.*f_I); %Depth averaged irradiance
            f_I2 = 1/H2*sum(dz.*f_I2);
            munet = Pmax*(f_I-r)/theta; %Depth averaged munet
            munet2 = Pmax2*(f_I2-r2)/theta2;

            term11 = (munet-ZP).*depaBio; %Growth Bio expression depth-averaged WC1 
            term12 = (munet2-ZP2).*depaBio2; %Growth Bio expression depth-averaged WC2 
            term21 = -BG*Biom(1,:); %BG Flux Bio expression depth-averaged WC1 
            term22 = -BG2*Biom2(1,:); %BG Flux Bio expression depth-averaged WC2
            term31 = mean(Ky)*ratio1*1/H*sum(dz.*diffBiom); %Exchange flux Bio expression depth-averaged WC1
            term32 = -mean(Ky)*ratio2*1/H*sum(dz.*diffBiom); %Exchange flux Bio expression depth-averaged WC2
            figure(18)
            sgtitle([filename, ' Run ',num2str(run_num)]) 
            subplot(2,2,1)
            plot(times(tstart:tend)/3600,term11(tstart:tend),times(tstart:tend)/3600,term12(tstart:tend))
            title('Depth-averaged Growth Flux')
            ylabel('(mg chl a/m^3)/s')
            xlabel('Time (hrs)')
            legend('WC 1','WC 2')
            subplot(2,2,2)
            plot(times(tstart:tend)/3600,term21(tstart:tend),times(tstart:tend)/3600,term22(tstart:tend))
            title('BG Flux')
            ylabel('(mg chl a/m^3)/s')
            xlabel('Time (hrs)')
            legend('WC 1','WC 2')
            subplot(2,2,3)
            plot(times(tstart:tend)/3600,term31(tstart:tend),times(tstart:tend)/3600,term32(tstart:tend))
            title('Depth-averaged Exchange Flux')
            ylabel('(mg chl a/m^3)/s') 
            xlabel('Time (hrs)')
            legend('WC 1','WC 2')
            subplot(2,2,4)
            plot(times(tstart:tend)/3600,term11(tstart:tend),times(tstart:tend)/3600,term12(tstart:tend))
            hold on
            plot(times(tstart:tend)/3600,term21(tstart:tend),times(tstart:tend)/3600,term22(tstart:tend))
            plot(times(tstart:tend)/3600,term31(tstart:tend),times(tstart:tend)/3600,term32(tstart:tend))
            hold off
            title('Depth Averaged Fluxes Bio')
            legend('WC 1 Growth','WC 2 Growth','WC 1 BG','WC 2 BG','WC 1 Exchange','WC 2 Exchange')
            ylabel('(mg chl a/m^3)/s') 
            xlabel('Time (hrs)')
        end

        if plotme(23) == 1 %Ky square waves
            figure(23)
            subplot(3,1,1)
            plot(times(tstart:tend)/3600/24,Kym(1,tstart:tend))
            title({[filename, ' Run ',num2str(run_num)]; 'Time Variable Ky'})
            ylabel('Effective Ky (m^2/s)')
            xlabel('Time (Days)')
            subplot(3,2,3)
            plot(times(tstart:tend)/3600/24,mean(Um(:,tstart:tend)))
            ylabel('Uave (m/s)')
            xlabel('Time (Days)')
            title('Spring-Neap WC1')
            subplot(3,2,5)
            plot(times(tstart:tend)/3600/24,mean(Um2(:,tstart:tend)))
            ylabel('Uave (m/s)')
            xlabel('Time (Days)')
            title('Spring-Neap WC2')
            subplot(3,2,4)
            plot(times/3600/24,mean(Pxm))
            xlabel('Time (Days)')
            ylabel('Ave Px (m/s^2)')
            title('Px')
            subplot(3,2,6)
            plot(times/3600/24,mean(Pxm2))
            xlabel('Time (Days)')
            ylabel('Ave Px2 (m/s^2)')
            title('Px2')
        end
        
        if plotme(33) == 1 %Both WCs Depth Averaged Salinity
%             figure(33)
%             subplot(1,2,1)
%             plot(times(tstart:tend)/3600/24,mean(Cm(:,tstart:tend)),'DisplayName','Tvar Ky')
%             ylabel('Salinity (psu)')
%             xlabel('Time (Days)')
%             title('Channel')
%             subplot(1,2,2)
%             plot(times(tstart:tend)/3600/24,mean(Cm2(:,tstart:tend)),'r-','DisplayName','Tvar Ky')
%             sgtitle({[filename, ' Run ',num2str(run_num)]; ['Depth Averaged Sal, Max Ky = ',num2str(Kytop), ' m^2/s']})
%             ylabel('Salinity (psu)')
%             xlabel('Time (Days)')
%             title('Shoal')
            
            figure(33)
            plot(times(teqstart:tend)/3600/24,(Cm(1,teqstart:tend)-Cm(N,teqstart:tend)),'DisplayName','Channel')
            hold on
           % plot(times(teqstart:tend)/3600/24,(Cm2(1,teqstart:tend)-Cm2(N2,teqstart:tend)),'r-','DisplayName','Shoal')
            hold off
           % title({[filename, ' Run ',num2str(run_num)]; ['Bottom-Top Salinity, Max Ky = ',num2str(Kytop), ' m^2/s']})
            ylabel('Channel Salinity Stratification (psu)')
            xlabel('Time (days)')
            %legend
        end
        
        if plotme(34) == 1 %Both WCs Depth Averaged Velocity
%             figure(34)
%             subplot(1,2,1)
%             plot(times(teqstart:tend)/3600/24,mean(Um(:,teqstart:tend)),'DisplayName','Tvar Ky')
%             ylabel('Velocity (m/s)')
%             xlabel('Time (Days)')
%             title('Channel')
%             subplot(1,2,2)
%             plot(times(teqstart:tend)/3600/24,mean(Um2(:,teqstart:tend)),'r-','DisplayName','Tvar Ky')
%             sgtitle({[filename, ' Run ',num2str(run_num)]; ['Depth Averaged Vel, Max Ky = ',num2str(Kytop), ' m^2/s']})
%             ylabel('Velocity (m/s)')
%             xlabel('Time (Days)')
%             title('Shoal')

            figure(34)
            plot(times(teqstart:tend)/3600/24,mean(Um(:,teqstart:tend)),'DisplayName','Channel')
            hold on
            plot(times(teqstart:tend)/3600/24,mean(Um2(:,teqstart:tend)),'r-','DisplayName','Shoal')
            hold off
        %    title({[filename, ' Run ',num2str(run_num)]; ['Depth Averaged Velocity, Max Ky = ',num2str(Kytop), ' m^2/s']})
            ylabel('Depth-Averaged Velocity (ms^{-1})')
            xlabel('Time (days)')
            legend
        end
        
        if plotme(36) == 1
            plot(times(tstart:tend)/3600/24,Biom2(N2,tstart:tend)-Biom(N,tstart:tend))
            title('Surface')
            ylabel('Biomass (mg chl a m^{-3})')
            legend('Shoal-Channel')
            xlabel('Time (days)')
        end
    end
end

%% Second group of plots
if whichplot == 2
    if plotme(17) == 1 %Add fit function and fill to depth averaged bio
        figure(17)
        subplot(1,2,1)
        hold on
        if phyvars(loadnum).tmax > tstart
            y = feval(c,times(tstart:phyvars(loadnum).tmax)/3600/24);
            plot(times(tstart:phyvars(loadnum).tmax)/3600/24,y(tstart:phyvars(loadnum).tmax),'b--','DisplayName','B1 Fit1')
            if wantfill == 1    
                ar = area(times(tstart:phyvars(loadnum).tmax)/3600/24,y); %FIXME - switch to just fill depabio - not the fit function?
                ar.FaceColor = 'blue';
                ar.DisplayName = 'Area B1 Fit1';
            end 
            if isfield(phyvars,'ff') == 0 && tend>phyvars(loadnum).tmax %FIXME: in future, determine in big loop when have these saved. For now have to calculate each loadum
                yy = fit(times(phyvars(loadnum).tmax:tend)',depaBio(phyvars(loadnum).tmax:tend)','exp1');
                test = coeffvalues(yy);
                test(2) = test(2)*3600*24;
                c = cfit(f,test(1),test(2));
                yy = feval(c,times(phyvars(loadnum).tmax:tend)/3600/24);
                plot(times(phyvars(loadnum).tmax:tend)/3600/24,yy,'b-.','DisplayName','B1 Fit2');
                if wantfill == 1
                    arr = area(times(phyvars(loadnum).tmax:tend)/3600/24,yy);
                    arr.FaceColor = 'blue';
                    arr.DisplayName = 'Area B1 Fit2';
                    disp(['Area B1 Approx: ', num2str(trapz(times(tstart:phyvars(loadnum).tmax)/3600/24,y)+trapz(times(phyvars(loadnum).tmax:tend)/3600/24,yy))])
                end 
            else
                if wantfill == 1
                    disp(['Area B1 Approx: ', num2str(trapz(times(tstart:phyvars(loadnum).tmax)/3600/24,y))])
                end
            end
        else
            y = feval(c,times/3600/24);
            plot(times/3600/24,y,'b--','DisplayName','B1 Fit')
            if wantfill == 1
                ar = area(times/3600/24,y);
                ar.FaceColor = 'blue';
                ar.DisplayName = 'Area B1';
                disp(['Area B1 Approx: ', num2str(trapz(times/3600/24,y))])
            end
        end
        legend
        hold off
        subplot(1,2,2)
        hold on
        if phyvars(loadnum).tmax2 > tstart 
            y2 = feval(c2,times(tstart:phyvars(loadnum).tmax2)/3600/24);
            plot(times(tstart:phyvars(loadnum).tmax2)/3600/24,y2(tstart:phyvars(loadnum).tmax2),'r--','DisplayName','B2 Fit1')
            if wantfill == 1
                ar2 = area(times(tstart:phyvars(loadnum).tmax2)/3600/24,y2);
                ar2.FaceColor = 'red';
                ar2.DisplayName = 'Area B2 Fit1';
            end
            if isfield(phyvars,'ff2') == 0 && tend>phyvars(loadnum).tmax2 
                yy2 = fit(times(phyvars(loadnum).tmax2:tend)',depaBio2(phyvars(loadnum).tmax2:tend)','exp1');
                test = coeffvalues(yy2);
                test(2) = test(2)*3600*24;
                c = cfit(f,test(1),test(2));
                yy2 = feval(c,times(phyvars(loadnum).tmax2:tend)/3600/24);
                plot(times(phyvars(loadnum).tmax2:tend)/3600/24,yy2,'r-.','DisplayName','B2 Fit2');
                if wantfill == 1 
                    arr2 = area(times(phyvars(loadnum).tmax2:tend)/3600/24,yy2);
                    arr2.FaceColor = 'red';
                    arr2.DisplayName = 'Area B2 Fit2';
                    disp(['Area B2 Approx: ', num2str(trapz(times(tstart:phyvars(loadnum).tmax2)/3600/24,y2)+trapz(times(phyvars(loadnum).tmax2:tend)/3600/24,yy2))])
                end
            else
                if wantfill == 1
                   disp(['Area B2 Approx: ', num2str(trapz(times(tstart:phyvars(loadnum).tmax2)/3600/24,y2))])
                end
            end
        else
            y2 = feval(c2,times/3600/24);
            plot(times/3600/24,y2,'r--','DisplayName','B2 Fit')
            if wantfill == 1 
                ar2 = area(times/3600/24,y2);
                ar2.FaceColor = 'red';
                ar2.DisplayName = 'Area B2';
                disp(['Area B2 Approx: ', num2str(trapz(times/3600/24,y2))])
            end
        end
        legend
        hold off
    end

    if plotme(24) == 1 && mod(run_num,BGtestend*BG2testend) == 0 && (BGtest > 1 && BG2test > 1 || mod(run_num,loadnum) == 0) %Can only plot when have all points and at least two points
        start = floor((run_num-1)/(BGtestend*BG2testend))*(BGtestend*BG2testend); %If starting from scratch
        figure(24)
        sgtitle([bigfile, ' Run ',num2str(run_num)]) 
        for i=1:loadnum-start
            BGx(i)=phyvars(start+i).condnum; %#ok<*SAGROW>
            BG2y(i)=phyvars(start+i).BG2test;
        end
        subplot(1,3,1)
        gs = gscatter(BGx,BG2y,full(bloom(start+1:loadnum)),[],'*',15);
        for i=1:length(gs)
            if gs(i).DisplayName == '0'
                gs(i).Color = 'r';
            elseif gs(i).DisplayName == '1'
                gs(i).Color = 'g';
            elseif gs(i).DisplayName == '3'
                gs(i).Color = 'b';
            end
        end
        xlabel('BG1 test #')
        ylabel('BG2 test #')
        title(['Bloom? WC1 Ky = ',num2str(Ky_range(Kytest))])
        subplot(1,3,2)
        gs = gscatter(BGx,BG2y,full(bloom2(start+1:loadnum)),[],'*',15);
         for i=1:length(gs)
            if gs(i).DisplayName == '0'
                gs(i).Color = 'r';
            elseif gs(i).DisplayName == '1'
                gs(i).Color = 'g';
            elseif gs(i).DisplayName == '7'
                gs(i).Color = 'y';
            end
        end
        xlabel('BG1 test #')
        ylabel('BG2 test #')
        title(['Bloom? WC2 Ky = ',num2str(Ky_range(Kytest))])
        subplot(1,3,3)
        gs=gscatter(BGx,BG2y,full(bothbloom(start+1:loadnum)),[],'*',15);
        for i=1:length(gs)
            if gs(i).DisplayName == '0'
                gs(i).Color = 'r';
                gs(i).DisplayName = 'Neither Bloom';
            elseif gs(i).DisplayName == '1'
                gs(i).Color = 'g';
                gs(i).DisplayName = 'WC2 bloom, WC1 only inc';
            elseif gs(i).DisplayName == '3'
                gs(i).Color = 'b';
                gs(i).DisplayName = 'WC1 only dec, WC2 only inc';
            elseif gs(i).DisplayName == '4'
                gs(i).Color = 'k';
                gs(i).DisplayName = 'Both only dec';
            elseif gs(i).DisplayName == '7'
                gs(i).Color = 'y';
                gs(i).DisplayName = 'WC2 dec only, WC1 only inc';
            elseif gs(i).DisplayName == '10'
                gs(i).Color = 'm';
                gs(i).DisplayName = 'Both bloom';
            elseif gs(i).DisplayName == '-2'
                gs(i).Color = 'c';
                gs(i).DisplayName = 'WC1 const dec WC2 bloom';
            end
        end
        xlabel('BG1 test #')
        ylabel('BG2 test #')
        title(['Bloom? WC2-WC1 Ky = ',num2str(Ky_range(Kytest))])
    else
        plotme(24) = 0; %Can't save if not have plot
    end

    if plotme(25) == 1 && mod(run_num,Kytestend*tidetestend) == 0 && (Kytest > 1 && condnum > 1 || mod(run_num,loadnum) == 0) %Can only plot when have FIXME - need all points? and at least two points
        start = floor((run_num-1)/(Kytestend*tidetestend))*(Kytestend*tidetestend); %If starting from scratch
        figure(25)
        sgtitle([bigfile, ' Run ',num2str(run_num)]) 
        for i=1:loadnum-start
            Kyx(i)=phyvars(start+i).Kytest;
            condy(i)=phyvars(start+i).condnum;
        end
        subplot(1,3,1)
        gs = gscatter(Kyx,condy,full(bloom(start+1:loadnum)),[],[],50);
        for i=1:length(gs)
            if gs(i).DisplayName == '0'
                gs(i).Color = 'r';
            elseif gs(i).DisplayName == '1'
                gs(i).Color = 'g';
            elseif gs(i).DisplayName == '3'
                gs(i).Color = 'b';
            end
        end
        xlabel('Ky test #')
        ylabel('BG2 test #')
        title(['Bloom? WC1 Ky = ',num2str(Ky_range(Kytest))])
        subplot(1,3,2)
        gs = gscatter(Kyx,condy,full(bloom2(start+1:loadnum)),[],[],50);
         for i=1:length(gs)
            if gs(i).DisplayName == '0'
                gs(i).Color = 'r';
            elseif gs(i).DisplayName == '1'
                gs(i).Color = 'g';
            elseif gs(i).DisplayName == '7'
                gs(i).Color = 'y';
            end
        end
        xlabel('Ky test #')
        ylabel('BG2 test #')
        title(['Bloom? WC2 Ky = ',num2str(Ky_range(Kytest))])
        subplot(1,3,3)
        gs=gscatter(Kyx,condy,full(bothbloom(start+1:loadnum)),[],[],50);
        for i=1:length(gs)
            if gs(i).DisplayName == '0'
                gs(i).Color = 'r';
                gs(i).DisplayName = 'Neither Bloom';
            elseif gs(i).DisplayName == '1'
                gs(i).Color = 'g';
                gs(i).DisplayName = 'WC2 bloom, WC1 only inc';
            elseif gs(i).DisplayName == '3'
                gs(i).Color = 'b';
                gs(i).DisplayName = 'WC1 only dec, WC2 only inc';
            elseif gs(i).DisplayName == '4'
                gs(i).Color = 'k';
                gs(i).DisplayName = 'Both only dec';
            elseif gs(i).DisplayName == '7'
                gs(i).Color = 'y';
                gs(i).DisplayName = 'WC2 dec only, WC1 only inc';
            elseif gs(i).DisplayName == '10' %#ok<*BDSCA>
                gs(i).Color = 'm';
                gs(i).DisplayName = 'Both bloom';
            elseif gs(i).DisplayName == '-2'
                gs(i).Color = 'c';
                gs(i).DisplayName = 'WC1 const dec WC2 bloom';
            elseif gs(i).DisplayName == '-1'
                gs(i).Color = [0.3010, 0.7450, 0.9330];
                gs(i).DisplayName = 'WC1 bloom only, WC2 only inc';
            elseif gs(i).DisplayName == '-3'
                gs(i).Color = [0.4940, 0.1840, 0.5560];
                gs(i).DisplayName = 'WC1 dec only, WC2 only inc';
            end
        end
        xlabel('Ky test #')
        ylabel('BG2 test #')
        title(['Bloom? WC2-WC1 Ky = ',num2str(Ky_range(Kytest))])
    else
        plotme(25) = 0; %Can't save if not have plot
    end

    if plotme(26) == 1 %Plot bar math terms WC1
        figure(26)
        if titletype == 1
            sgtitle(['Cond 2 Math Terms WC1 Run ', num2str(rang(:,1)'),' (Kymax = ',num2str(vars(rang(1,1)).Kytop),')']) 
        else
            sgtitle(['Cond 2 Math Terms WC1 Run ', num2str(rang(:,1)'),' (Days = ',num2str(vars(rang(1,1)).numdays),')']) 
        end
        subplot(3,2,1)
        b = bar(X,bmunetbB(rang));
        b(2).FaceColor='k';
        legend({'Tvar','Constant'},'Location','northwestoutside')
        title('$\overline{\mu_{net}}\overline{B}$','Interpreter','Latex')
        subplot(3,2,2)
        b = bar(X,bmunetprBpr(rang));
        b(2).FaceColor='k';
        title('$\overline{\mu_{net}''B''}$','Interpreter','Latex')
        subplot(3,2,3)
        b = bar(X,bKybdB(rang));
        b(2).FaceColor='k';
        title('$\overline{K_y}\overline{dB}$','Interpreter','Latex')
        subplot(3,2,4)
        b = bar(X,bKyprdBpr(rang));
        b(2).FaceColor='k';
        title('$\overline{K_y''dB''}$','Interpreter','Latex')
        subplot(3,2,5)
        b = bar(X,bioexp(rang));
        b(2).FaceColor='k';
        title('$-(ZP+\frac{BG}{H})\overline{B}$','Interpreter','Latex')
        subplot(3,2,6)
        b = bar(X,total(rang));
        b(2).FaceColor='k';
        title('$\frac{\partial\overline{B}}{\partial t}$','Interpreter','Latex')
    else
        plotme(26) = 0; %Can't save if not have plot
    end
    
    if plotme(27) == 1 %Plot bar math terms WC2
        figure(27)
        if titletype == 1
            sgtitle(['Cond 2 Math Terms WC1 Run ', num2str(rang(:,1)'),' (Kymax = ',num2str(vars(rang(1,1)).Kytop),')']) 
        else
            sgtitle(['Cond 2 Math Terms WC1 Run ', num2str(rang(:,1)'),' (Days = ',num2str(vars(rang(1,1)).numdays),')']) 
        end
        subplot(3,2,1)
        b=bar(X,bmunet2bB2(rang));
        b(1).FaceColor = 'r';
        b(2).FaceColor='k';
        legend({'Tvar','Constant'},'Location','northwestoutside')
        title('$\overline{\mu_{net2}}\overline{B_2}$','Interpreter','Latex')
        subplot(3,2,2)
        b=bar(X,bmunet2prB2pr(rang));
        b(1).FaceColor = 'r';
        b(2).FaceColor='k';
        title('$\overline{\mu_{net2}''B_2''}$','Interpreter','Latex')
        subplot(3,2,3)
        b=bar(X,bKybdB2(rang));
        b(1).FaceColor = 'r';
        b(2).FaceColor='k';
        title('$\overline{K_y}\overline{dB}$','Interpreter','Latex')
        subplot(3,2,4)
        b=bar(X,bKyprdBpr2(rang));
        b(1).FaceColor = 'r';
        b(2).FaceColor='k';
        title('$\overline{K_y''dB''}$','Interpreter','Latex')
        subplot(3,2,5)
        b=bar(X,bioexp2(rang));
        b(1).FaceColor = 'r';
        b(2).FaceColor='k';
        title('$-(ZP_2+\frac{BG_2}{H_2})\overline{B_2}$','Interpreter','Latex')
        subplot(3,2,6)
        b=bar(X,total2(rang));
        b(1).FaceColor = 'r';
        b(2).FaceColor='k';
        title('$\frac{\partial\overline{B_2}}{\partial t}$','Interpreter','Latex')
    else
        plotme(27) = 0; %Can't save if not have plot
    end
    
    if plotme(28) == 1 %Plot math terms scatter plot WC1
        figure(28)
        sgtitle('Cond 2 Math Terms WC1')
        [x,~,y]=find(bmunetbB);
        [~,~,c]=find(numdays(x));
        subplot(3,2,1)
        gs = gscatter(x,y,c);
%         for i=1:length(gs) %FIXME: different marker type even vs odd?
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        title('$\overline{\mu_{net}}\overline{B}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        subplot(3,2,2)
        [x,~,y]=find(bmunetprBpr);
        [~,~,c]=find(numdays(x));
        gs = gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        title('$\overline{\mu_{net}''B''}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
        subplot(3,2,3)
        [x,~,y]=find(bKybdB);
        [~,~,c]=find(numdays(x));
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        title('$\overline{K_y}\overline{dB}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
        subplot(3,2,4)
        [x,~,y]=find(bKyprdBpr);
        [~,~,c]=find(numdays(x));
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        title('$\overline{K_y''dB''}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
        subplot(3,2,5)
        [x,~,y]=find(bioexp);
        [~,~,c]=find(numdays(x));
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        title('$-(ZP+\frac{BG}{H})\overline{B}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
        subplot(3,2,6)
        [x,~,y]=find(total);
        [~,~,c]=find(numdays(x));
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        legend off
        title('$\frac{\partial\overline{B}}{\partial t}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
    else
        plotme(28) = 0; %Can't save if not have plot
    end
    
    if plotme(29) == 1 %Plot math terms scatter plot WC2
        figure(29)
        sgtitle('Cond 2 Math Terms WC2')
        [x,~,y]=find(bmunet2bB2);
        [~,~,c]=find(numdays(x));
        subplot(3,2,1)
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        title('$\overline{\mu_{net2}}\overline{B_2}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        subplot(3,2,2)
        [x,~,y]=find(bmunet2prB2pr);
        [~,~,c]=find(numdays(x));
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        title('$\overline{\mu_{net2}''B_2''}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
        subplot(3,2,3)
        [x,~,y]=find(bKybdB2);
        [~,~,c]=find(numdays(x));
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        title('$\overline{K_y}\overline{dB}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
        subplot(3,2,4)
        [x,~,y]=find(bKyprdBpr2);
        [~,~,c]=find(numdays(x));
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        title('$\overline{K_y''dB''}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
        subplot(3,2,5)
        [x,~,y]=find(bioexp2);
        [~,~,c]=find(numdays(x));
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        title('$-(ZP_2+\frac{BG_2}{H_2})\overline{B_2}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
        subplot(3,2,6)
        [x,~,y]=find(total2);
        [~,~,c]=find(numdays(x));
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '3'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             end
%         end
        title('$\frac{\partial\overline{B_2}}{\partial t}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
    else
        plotme(29) = 0; %Can't save if not have plot
    end
    
    if plotme(30) == 1 %Plot scatter plot only barKy'dB' terms FIXME - see live script. Change way plot?
        figure(30)
        subplot(3,2,1)
        x = 11:2:90;
        y = bKyprdBpr(x);
        [~,~,c]=find(numdays(x));
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '28'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             elseif gs(i).DisplayName == '56'
%                 gs(i).Color = [0 0 1];
%             end
%         end
        hold on
        for i=1:Kytestend-1
            x1 = x(1+BGtestend*(i-1):BGtestend+BGtestend*(i-1));
            model = polyfit(x1,y(1+BGtestend*(i-1):BGtestend+BGtestend*(i-1)),1);
            slope(i) = model(1);
            y1 = polyval(model,x1);
            plot(x1,y1,'r-','HandleVisibility','off');
            txt = ['Slope = ' num2str(slope(i))];
         %   text(x1(1),y1(1)-0.3*10^-4,txt);
        end
        hold off
        title('Channel $\overline{K_y''dB''}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        
        subplot(3,2,3)
        x = Ky_range(2:Kytestend);
        y = slope;
        plot(x,y,'o');
        title('Channel','Interpreter','Latex')
        hold on
     %   model = polyfit(x,y,1);
   %     y1 = polyval(model,x);
    %    plot(x,y1,'r-','HandleVisibility','off');
        txt = ['Slope = ' num2str(model(1))];
     %   text(x(1),y1(1),txt);
        xlabel('K_y max')
        ylabel('Slope')
        hold off
        
        subplot(3,2,2)
        x = 11:2:90;
        y = bKyprdBpr2(x);
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '28'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             elseif gs(i).DisplayName == '56'
%                 gs(i).Color = [0 0 1];
%             end
%         end
        hold on
        for i=1:Kytestend-1
            x1 = x(1+BGtestend*(i-1):BGtestend+BGtestend*(i-1));
            model = polyfit(x1,y(1+BGtestend*(i-1):BGtestend+BGtestend*(i-1)),1);
            slope2(i) = model(1);
            y1 = polyval(model,x1);
            plot(x1,y1,'r-');
            txt = ['Slope = ' num2str(slope2(i))];
         %   text(x1(1),y1(1)-0.1*10^-4,txt);
        end
        hold off
        title('Shoal $\overline{K_y''dB''}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
        
        subplot(3,2,4)
        x = Ky_range(2:Kytestend);
        y = slope2;
        plot(x,y,'ro');
        title('Shoal','Interpreter','Latex')
        hold on
%         model = polyfit(x,y,1);
%         y1 = polyval(model,x);
 %       plot(x,y1,'r-','HandleVisibility','off');
        txt = ['Slope = ' num2str(model(1))];
      %  text(x(1),y1(1),txt);
        xlabel('K_y max')
        ylabel('Slope')
        hold off
        
        subplot(3,2,5)
        x = 11:2:90;
        y = bKyprdBpr(x)./bKybdB(x);
        [~,~,c]=find(numdays(x));
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '28'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             elseif gs(i).DisplayName == '56'
%                 gs(i).Color = [0 0 1];
%             end
%         end
        title('Channel $\frac{\overline{K_y''dB''}}{\overline{K_ydB}}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
        
        subplot(3,2,6)
        x = 11:2:90;
        y = bKyprdBpr2(x)./bKybdB2(x);
        gs=gscatter(x,y,c);
%         for i=1:length(gs)
%             if gs(i).DisplayName == '7'
%                 gs(i).Color = [0.4940 0.1840 0.5560];
%             elseif gs(i).DisplayName == '14'
%                 gs(i).Color = [0.4660 0.6740 0.1880];
%             elseif gs(i).DisplayName == '28'
%                 gs(i).Color = [0.9290 0.6940 0.1250];
%             elseif gs(i).DisplayName == '56'
%                 gs(i).Color = [0 0 1];
%             end
%         end
        title('Shoal $\frac{\overline{K_y''dB''}}{\overline{K_ydB}}$','Interpreter','Latex')
        xlabel('Run #')
        ylabel(' ')
        legend off
    end
    
    if plotme(31) == 1 %Plot bar math terms relative
        figure(31)
        if titletype == 1
            sgtitle(['Cond 2 Math Terms (Kymax = ',num2str(vars(rang(1,1)).Kytop),')']) %Run ', num2str(rang(:,1)'),' 
        else
            sgtitle(['Cond 2 Math Terms Run ', num2str(rang(:,1)'),' (Days = ',num2str(vars(rang(1,1)).numdays),')']) 
        end
        subplot(3,3,1)
        b = bar(X,(bmunetbB(rang)+bmunetprBpr(rang))./total(rang));
        b(2).FaceColor='k';
     %   legend({'Tvar','Constant'},'Location','northwestoutside')
        title('$\frac{\overline{\mu_{net}}\overline{B}+\overline{\mu_{net}''B''}}{\frac{\partial\overline{B}}{\partial t}}$','Interpreter','Latex')
        subplot(3,3,4)
        b = bar(X,bmunetprBpr(rang)./bmunetbB(rang));
        b(2).FaceColor='k';
        title('$\frac{\overline{\mu_{net}''B''}}{\overline{\mu_{net}}\overline{B}}$','Interpreter','Latex')
        subplot(3,3,2)
        b=bar(X,(bmunet2bB2(rang)+bmunet2prB2pr(rang))./total2(rang));
        b(1).FaceColor = 'r';
        b(2).FaceColor='k';
   %     legend({'Tvar','Constant'},'Location','northwestoutside')
        title('$\frac{\overline{\mu_{net2}}\overline{B_2}+\overline{\mu_{net2}''B_2''}}{\frac{\partial\overline{B_2}}{\partial t}}$','Interpreter','Latex')
        subplot(3,3,5)
        b=bar(X,bmunet2prB2pr(rang)./bmunet2bB2(rang));
        b(1).FaceColor = 'r';
        b(2).FaceColor='k';
        title('$\frac{\overline{\mu_{net2}''B_2''}}{\overline{\mu_{net2}}\overline{B_2}}$','Interpreter','Latex')
        subplot(3,3,9)
        b = bar(X,bKyprdBpr(rang)./bKybdB(rang));
        b(1).FaceColor = [0.5 0 0.5];
        b(2).FaceColor='k';
     %   legend('Both WCs','Location','northwestoutside')
        title('$\frac{\overline{K_y''dB''}}{\overline{K_y}\overline{dB}}=R_k$','Interpreter','Latex')
        subplot(3,3,8)
        b=bar(X,(bKybdB2(rang)+bKyprdBpr2(rang))./total2(rang));
        b(1).FaceColor = 'r';
        b(2).FaceColor='k';
        title('$\frac{\overline{K_y}\overline{dB}+\overline{K_y''dB''}}{\frac{\partial\overline{B_2}}{\partial t}}$','Interpreter','Latex')
        subplot(3,3,6)
        b=bar(X,bioexp2(rang)./total2(rang));
        b(1).FaceColor = 'r';
        b(2).FaceColor='k';
        title('$\frac{-(ZP_2+\frac{BG_2}{H_2})\overline{B_2}}{\frac{\partial\overline{B_2}}{\partial t}}$','Interpreter','Latex')
        subplot(3,3,3)
        b = bar(X,bioexp(rang)./total(rang));
        b(2).FaceColor='k';
        title('$\frac{-(ZP+\frac{BG}{H})\overline{B}}{\frac{\partial\overline{B}}{\partial t}}$','Interpreter','Latex')
        subplot(3,3,7)
        b=bar(X,(bKybdB(rang)+bKyprdBpr(rang))./total(rang));
        b(2).FaceColor='k';
        title('$\frac{\overline{K_y}\overline{dB}+\overline{K_y''dB''}}{\frac{\partial\overline{B}}{\partial t}}$','Interpreter','Latex')
    else
        plotme(31) = 0; %Can't save if not have plot
    end
    
    if plotme(32) == 1 %Plot contour plots against Kymax, Numdays 
        figure(32)
        subplot(1,2,1) %Kyave = 5
        x = 11:2:50;
        bottom = min(min(bKyprdBpr./bKybdB));
        top = max(max(bKyprdBpr./bKybdB));
        for i=1:length(x)
            numx(i) = num_range(phyvars(x(i)).BGtest)*2;
            Kyy(i) = Ky_range(phyvars(x(i)).Kytest);
            Z = bKyprdBpr(x)./bKybdB(x);
        end
        X = linspace(numx(1),numx(BGtestend),100);
        Y = linspace(Kyy(1),Kyy((Kytestend-1)/2*BGtestend),80);
        [numxq, Kyyq] = meshgrid(X,Y);
        numxq = numxq';
        Kyyq = Kyyq';
        fitobject = fit([numx',Kyy'],Z,'poly21');
        contourf(reshape(numx,5,4),reshape(Kyy,5,4),reshape(Z,5,4),'LineStyle','none');
        hold on
        contour(numxq,Kyyq,fitobject(numxq,Kyyq),'k-')
        hold off
        if wantsamecolor == 1 
            caxis manual
            caxis([bottom top]);
        end
        xlabel('Period (days)')
        ylabel('Kymax (m^2/s)')
        title('Kyave = 5 m^2/s')
        h = colorbar;
        xlabel(h,'R_k')
        colormap(colorscheme)
        subplot(1,2,2) %Kyave = 25
        x = 51:2:90;
        for i=1:length(x)
            numx(i) = num_range(phyvars(x(i)).BGtest)*2;
            Kyy(i) = Ky_range(phyvars(x(i)).Kytest);
            Z = bKyprdBpr(x)./bKybdB(x);
        end
        X = linspace(numx(1),numx(BGtestend),100);
        Y = linspace(Kyy(1),Kyy((Kytestend-1)/2*BGtestend),80);
        [numxq, Kyyq] = meshgrid(X,Y);
        numxq = numxq';
        Kyyq = Kyyq';
        fitobject2 = fit([numx',Kyy'],Z,'poly21');
        contourf(reshape(numx,5,4),reshape(Kyy,5,4),reshape(Z,5,4),'LineStyle','none');
        shading interp
        if wantsamecolor == 1 
            caxis manual
            caxis([bottom top]);
        end
        colormap(colorscheme)
        xlabel('Period (days)')
        ylabel('Kymax (m^2/s)')
        title('Kyave = 25 m^2/s')
        h = colorbar;
        xlabel(h,'R_k')
        hold on
        contour(numxq,Kyyq,fitobject2(numxq,Kyyq),'k-')
        hold off
    else
        plotme(32) = 0; %Can't save if not have plot
    end
end
