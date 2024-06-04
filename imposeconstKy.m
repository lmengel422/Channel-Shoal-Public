% ***************************************************************************
%  Code to make plot 35; impose black line of a different depth-averaged
%   bio with time constant Ky against the time variable version in plot 17. 
% ***************************************************************************

figure(17)
savefig %Temporary file of previous biology timeseries 

%Run code with new constant Ky value
if loadmePhy == 1 
    Kyave = phyvars(run_num).Kyave; %Take from loaded or recently saved file
elseif saveme == 1
    Kyave = phyvars(run_num).Kyave; %Take from loaded or recently saved file
    loadnum = run_num; %Use first run's value to create plots, but don't affect newly imposed run
else
    loadnum = run_num; %Use first run's value to create plots, but don't affect newly imposed run
end
clearvars -except iii Kytest Kytestend ii loadme loadnum tidetestend jj loadmePhy kk daystest daystestend ll eftest eftestend condnum Kyave imposeKy 
watercolumn
calc_phyv2  

close all
openfig('Untitled.fig'); %Bring up temp fig

%% Plot on saved fig
figure(1)
subplot(1,2,1)
hold on
plot(times(tstart:tend)/3600/24,depaBio(tstart:tend),'k-','DisplayName',['Const Ky = ', num2str(Kyave),' m^2/s']) %Add imposed Ky bio trajectory
if teqstart > 1 || condnum == 2 %Plot average of certain section. Only enabled when zoom in currently. 
    plot(times(teqstart:tend)/3600/24,phyvars(loadnum).Bbar*ones(length(teqstart:tend),1),'DisplayName','Tvar')
    txt = ['$\overline{B_1}$ = ' num2str(phyvars(loadnum).Bbar)];
    text(times(teqstart)/3600/24,phyvars(loadnum).Bbar-1,txt,'Interpreter','Latex');
    plot(times(teqstart:tend)/3600/24,mean(depaBio(teqstart:tend))*ones(length(teqstart:tend),1),'DisplayName','Const')
    txt = ['$\overline{B_1}$ = ' num2str(mean(depaBio(teqstart:tend)))];
    text(times(teqstart)/3600/24,mean(depaBio(teqstart:tend))+1,txt,'Interpreter','Latex');
end
% legend('Location','south')
subplot(1,2,2)
hold on
plot(times(tstart:tend)/3600/24,depaBio2(tstart:tend),'k-','DisplayName',['Const Ky = ', num2str(Kyave),' m^2/s'])
if teqstart > 1 || condnum == 2
    plot(times(teqstart:tend)/3600/24,phyvars(loadnum).Bbar2*ones(length(teqstart:tend),1),'DisplayName','Tvar')
    txt = ['$\overline{B_2}$ = ' num2str(phyvars(loadnum).Bbar2)];
    text(times(teqstart)/3600/24,phyvars(loadnum).Bbar2+3,txt,'Interpreter','Latex');
    plot(times(teqstart:tend)/3600/24,mean(depaBio2(teqstart:tend))*ones(length(teqstart:tend),1),'DisplayName','Const')
    txt = ['$\overline{B_2}$ = ' num2str(mean(depaBio2(teqstart:tend)))];
    text(times(teqstart)/3600/24,mean(depaBio2(teqstart:tend))-3,txt,'Interpreter','Latex');
end
% legend('Location','south')

%% Save this figure
if (saveme == 1 && loadmePhy == 0) || savefigme == 1
   saveas(figure(1),[savefile,'\',filename, ' Run ',num2str(run_num-1),' Plot ',num2str(35),'.png'])
   disp(['Plot 35 saved. Run_num: ',num2str(run_num)])
end

delete('Untitled.fig') %Get rid of temporary file