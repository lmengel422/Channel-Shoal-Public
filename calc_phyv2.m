% *******************************************************************
%  Output values to calculate/interpret. v2: includes mutilBtil.
% *******************************************************************
if loadmePhy ~= 1
    if first_water_column == 1
        [maxBio,tmaxmax] = max(max(Biom(:,tstart:tend))); %Maximum biomass WC1(represent bloom) and where it happens
        dtimavKz = mean(mean(Kzm(:,tstart:tend))); %Depth and time averaged Kz WC1
        [maxdepBio,tmax] = max(depaBio); %Max depth integrated B WC1
        if tmax > tstart && tmax<tend
            f = fit(times(tstart:tmax)',depaBio(tstart:tmax)','exp1'); %Exponential growth fit to depth avg max WC1. FIXME: Matlab 2020b no longer allows Curve fitting toolbox for free.
            ff = fit(times(tmax:tend)',depaBio(tmax:tend)','exp1'); %Exponential decay to end WC1
        else
            f = fit(times(tstart:tend)',depaBio(tstart:tend)','exp1'); %Exponential decay to end WC1
            ff = 0;
        end
        %Equilibrium calculations
        Bbar = mean(depaBio(teqstart:tend));
        munetave = mean(mean(munetm(:,teqstart:tend)));
        munetB = munetave * Bbar;
        muprBpr = mean(mean(munetm(:,teqstart:tend)-munetave).*mean(Biom(:,teqstart:tend)-Bbar)); %Mean averages each column of a matrix (depth), then total of array (time)
        mubar = mean(munetm(:,teqstart:tend)); %Depth-average munet equilibrium
        mutilBtil = mean((munetm(:,teqstart:tend)-mubar).*(Biom(:,teqstart:tend)-depaBio(teqstart:tend))); %mutilBtil
        B0bar = mean(Biom(1,teqstart:tend)); %time-average bottom B
    end
    if second_water_column == 1
        [maxBio2,tmaxmax2] = max(max(Biom2(:,tstart:tend))); %Max biomass WC2 and where it happens
        dtimavKz2 = mean(mean(Kzm2(:,tstart:tend))); %Depth and time averaged Kz WC2
        [maxdepBio2,tmax2] = max(depaBio2); %Max depth integrated B WC2
        if tmax2 > tstart && tmax2 < tend
            f2 = fit(times(tstart:tmax2)',depaBio2(tstart:tmax2)','exp1'); %Exponential growth fit up to depth avg max WC2 FIXME: Matlab 2020b no longer allows Curve fitting toolbox for free.
            ff2 = fit(times(tmax2:tend)',depaBio(tmax2:tend)','exp1'); %Exponential decay to end WC2
        else
            f2 = fit(times(tstart:tend)',depaBio2(tstart:tend)','exp1'); %Exponential decay fit up to end WC2
            ff2 = 0;
        end
        %Equilibrium calculations
        Bbar2 = mean(depaBio2(teqstart:tend));
        munetave2 = mean(mean(munetm2(:,teqstart:tend)));
        munetB2 = munetave2 * Bbar2;
        muprBpr2 = mean(mean(munetm2(:,teqstart:tend)-munetave2).*mean(Biom2(:,teqstart:tend)-Bbar2)); %Mean averages each column if a matrix (depth), then total if array (time)
        mubar2 = mean(munetm2(:,teqstart:tend)); %depth-average mutnet equilibrium
        mutilBtil2 = mean((munetm2(:,teqstart:tend)-mubar2).*(Biom2(:,teqstart:tend)-depaBio2(teqstart:tend))); %mutilBtil
        B0bar2 = mean(Biom2(1,teqstart:tend)); %time-average bottom B
    end
    if first_water_column == 1 && second_water_column == 1 
        %Equilibrium calculations
        Kyave = mean(mean(Kym(:,teqstart:tend)));
        dBiomave = mean(mean(diffBiom(:,teqstart:tend))); %Only difference between where the two water columns exchange. 
        KydB = Kyave * dBiomave;
        KyprdBpr = mean(mean(Kym(:,teqstart:tend)-Kyave).*mean(diffBiom(:,teqstart:tend)-dBiomave)); %Mean averages each column if a matrix (depth), then total if array (time)
    end
end

bigfile = ['PhytoBloom ',filename]; %Where to save the "post" data. Related to filename
savewhat = {'run_num','Kytest','condnum','daystest','eftest','maxBio','tmaxmax','dtimavKz',...
    'maxdepBio','tmax','f','ff','maxBio2','tmaxmax2','dtimavKz2',...
    'maxdepBio2','tmax2','f2','ff2','teqstart','Bbar','munetave','munetB',...
    'muprBpr','Bbar2','munetave2','munetB2','muprBpr2','Kyave','dBiomave',...
    'KydB','KyprdBpr','mutilBtil','mutilBtil2','B0bar','B0bar2'}; %Which of the calculated variables want to save

% *******************************************************************
%  Save calculated values to .MAT file 
% *******************************************************************
if saveme == 1 && loadmePhy == 0 %Whether or not to save this run
    if isfile([savefile,'\',bigfile,'.mat']) %See if this file exists
        mmm = matfile([savefile,'\',bigfile,'.mat']); %Evaluate details of file without loading: save time

        if isprop(mmm,'phyvars')==1  %By third run, phyvars exist. Can just load 
            load([savefile,'\',bigfile,'.mat']); 
        else
            phyvars = load([savefile,'\',bigfile,'.mat']); %Create struct of phyvars
        end

        if allvars == 1
            save('lastmax')      %Save initial variables this run
        else
            for i = 1:numel(savewhat)
                if i == 1
                    save('lastmax', savewhat{i}); %Write file to start
                else
                    save('lastmax', savewhat{i},'-append'); %Save variables saved in all other tests
                end
            end  
        end

        thisrun = load('lastmax');

        phyvars = [phyvars,thisrun]; %#ok<*AGROW>
        save([savefile,'\',bigfile,'.mat'],'phyvars'); %Save new variable set

        if export == 1 %Export struct to excel
            writetable(struct2table(phyvars),[savefile,'\',bigfile,'.xlsx']); 
        end
        delete('lastmax.mat')
    else %First run save new file
        if allvars == 1
            save([savefile,'\',bigfile,'.mat'])      %Save initial variables this run
        else
            for i = 1:numel(savewhat)
                if i == 1
                    save([savefile,'\',bigfile,'.mat'], savewhat{i}); %Write file to start
                else
                    save([savefile,'\',bigfile,'.mat'], savewhat{i},'-append'); %Save variables saved in all other tests
                end
            end  
        end
        phyvars = load([savefile,'\',bigfile,'.mat']);
    end
    disp(['File saved. File: ', bigfile, '. Run_num: ' num2str(run_num)])
end

% *******************************************************************
%  Load calculated values from a .MAT file (replace what just made).
% *******************************************************************
if loadmePhy == 1 && saveme == 0 %Whether or not to take calculated variables from loaded file. Cannot be loaded if saved same run.
    mmm = matfile([savefile,'\',bigfile,'.mat']); %Evaluate details of file without loading: save time

    if isprop(mmm,'phyvars')==1  %By third run, phyvars exist. Can just load 
        load([savefile,'\',bigfile,'.mat']); 
    else
        phyvars = load([savefile,'\',bigfile,'.mat']); %Create struct of phyvars
    end
    run_num = loadnum; %Make easy for following plots

    names = fieldnames(phyvars); %Overwrite force variables for plot titles
    for i = 1:4 %Only changing 4 variables (Ky, Tide, days, ef)
        assignin('caller', names{i}, phyvars(loadnum).(names{i}));
    end
    disp(['File loaded. File: ', bigfile, '. Loadnum: ',num2str(loadnum)])
end

% *******************************************************************
%  Save desired plots to .fig  or .png file FIXME - move?
% *******************************************************************
if (savefigme == 1 && loadmePhy == 0)
    if isfolder(filename) == 0 %Make directory to save plot outputs
        mkdir(filename) 
    end
    j=find(plotme==1); %only save plots that were desired/made each run
    if isempty(j) ~= 1 %If there are any figures to save
        for i=1:length(j)
            if wantfig == 1 && wantpng == 0
                h(i) = figure(j(i)); %#ok<SAGROW> %Make array of figures 
            elseif wantpng == 1 && wantfig == 0 %Saveas png
                saveas(figure(j(i)),[savefile,'\',filename, ' Run ',num2str(run_num),' Plot ',num2str(j(i)),'.png']) %Save from this run to one fig file. Have to compbine later
            end
        end
        if wantfig == 1 && wantpng == 0 %Save as fig
            savefig(h,[savefile,'\',filename, ' Run ',num2str(run_num),'.fig'],'compact') %Save from this run to one fig file. Have to compbine later
        end
    end
    disp(['Figures ',num2str(j),' saved. Run_num: ',num2str(run_num)])
end