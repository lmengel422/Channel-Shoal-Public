% ********************************************************************
% Runs watercolumn and interprets/saves important data outputs from it.
%       If want to start from watercolumn instead, need to provide values
%       for force variables directly below but before loop.
%       Total other files needed: watercolumn.m, wc_advance.m, wc_setup.m,
%       wc_preallocate.m, plots.m, calc_phyv2.m
%       Fall 2019 & beyond: Lily Engel
% ********************************************************************

% *******************************************************************
%%  Variables to force code (force variables)
% *******************************************************************
Kytest = 0; %Starting value for Ky test. Should be 0 so can start at test=1 in loop
daystest = 1; %Starting value for numdays test. Should be 1 so can start at daystest=2 in loop (1st value is dummy variable)
eftest = 1; %Starting value for ebbfl test. Should be 1 so can start at eftest=2 in loop
Kytestend = 9; %Length Ky vector going to loop through 9
tidetestend = 4; %Number Ky conditions going to loop through 4
daystestend = 9; %Length numdays vector going to loop through. One less than total because of dummy.  9
eftestend = 2; %Length ebb/flood  vector going to loop through. One less than total because of dummy.
loadme = 1; %Whether or not to recreate a condition from a previous test
loadmePhy = 1; %Whether or not to recreate phytobloom test; don't run watercolumn again
loadnum = 273; %Which row of the loaded struct to recreate. Should start at 1 if looping through watercolumn or final number investigating if loadmePhy. 269
condnum = 0; %Which timevariable Ky expression to use. Start at 0 if going to loop.
    %Set to 1 for tidally dependent square wave
    %Set to 2 to have on numdays days, off next
    %Set to 3 to determine based on some U value
    %Set to 4 to have condition 1 and 3 both at same time
imposeKy = 0; %Start not imposing Ky, then after first loop do so
Kyave = 0; %Variable needed for second run later

% *******************************************************************
%%  Loop through test variables; run watercolumn for each test case
% *******************************************************************
for ii = 1:Kytestend
    clearvars -except Kytest Kytestend ii loadme loadnum tidetestend jj loadmePhy kk daystest daystestend ll eftest eftestend condnum Kyave imposeKy %Clear all except force variables
    Kytest = Kytest + 1; %Which test this is
    
    for jj = 1:tidetestend
        condnum = condnum + 1;
        
        if condnum == 1 || condnum == 4
            for ll = 1:eftestend
                tic
                eftest = eftest + 1; %Using this variable for ebb or flood FIXME
                % *******************************************************************
                %  Call code to run model with current test scenario
                % *******************************************************************
                watercolumn
                calc_phyv2    
                if plotme(35) == 2 && (mod(loadnum,2) == 1 || loadme == 0) %Only impose if requested and correct line in saved file
                    imposeKy = 1; %#ok<NASGU>
                    imposeconstKy
                    imposeKy = 0;
                end
                toc

                if loadmePhy == 1 || loadme == 1 %Don't need to go through loops
                    break
                end
            end
            eftest = 1; %Reset loop
            daystest = 1;
        elseif condnum == 2
            for kk = 1:daystestend
                tic
                daystest = daystest + 1; %Using this variable for number of days 
                if daystest == 2 || daystest == 3 || daystest == 4 || daystest == 5 %If numdays = 7 (period is 14 days)
                    for ll = 1:eftestend
                       eftest = eftest + 1; %Use ebb/flood variable to determine if in or out of sync with spring neap
                       watercolumn
                       calc_phyv2
                        if plotme(35) == 2 && (mod(loadnum,2) == 1 || loadme == 0)
                            imposeKy = 1; %#ok<NASGU>
                            imposeconstKy
                            imposeKy = 0;
                        end
                    end
                eftest = 1; %reset loop
                else
                    watercolumn
                    calc_phyv2
                    if plotme(35) == 2 && (mod(loadnum,2) == 1 || loadme == 0)
                        imposeKy = 1; %#ok<NASGU>
                        imposeconstKy
                        imposeKy = 0;
                    end
                end
                toc

                if loadmePhy == 1 || loadme == 1 %Don't need to go through loops
                    break
                end
            end
            daystest = 1; %Reset loop
            eftest = 1;
        else
            tic
            watercolumn
            calc_phyv2
            if plotme(35) == 2 && (mod(loadnum,2) == 1 || loadme == 0)
                imposeKy = 1; %#ok<NASGU>
                imposeconstKy
                imposeKy = 0;
            end
            daystest = 1; %Reset loop
            eftest = 1;
            toc
        end
            
        if loadmePhy == 1 || loadme == 1 %Don't need to go through loops
            break
        end
    end
    condnum = 0; %Reset loop
    if loadmePhy == 1 || loadme == 1 %Don't need to go through loops
        break
    end
end