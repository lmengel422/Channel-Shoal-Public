% *******************************************************************
%%  Code to interpret results from struct "phyvars"
% *******************************************************************

bmunetbB = sparse(length(phyvars),1);
bmunetprBpr = sparse(length(phyvars),1);
bmunet2bB2 = sparse(length(phyvars),1);
bmunet2prB2pr = sparse(length(phyvars),1);
bKybdB = sparse(length(phyvars),1);
bKyprdBpr = sparse(length(phyvars),1);
bKybdB2 = sparse(length(phyvars),1);
bKyprdBpr2 = sparse(length(phyvars),1);
bioexp = sparse(length(phyvars),1);
bioexp2 = sparse(length(phyvars),1);
total = sparse(length(phyvars),1);
total2 = sparse(length(phyvars),1);
numdays = sparse(length(vars),1);
Kymax = sparse(length(vars),1);
bmutilBtil = sparse(length(phyvars),1);
bmutilBtil2 = sparse(length(phyvars),1);
Rk = sparse(length(phyvars),1);
Bbar = sparse(length(phyvars),1);
Bbar2 = sparse(length(phyvars),1);
B0bar = sparse(length(phyvars),1);
B0bar2 = sparse(length(phyvars),1);

for i=1:length(phyvars)
    bmunetbB(i) = phyvars(i).munetB;
    bmunetprBpr(i) = phyvars(i).muprBpr;
    bmunet2bB2(i) = phyvars(i).munetB2;
    bmunet2prB2pr(i) = phyvars(i).muprBpr2;
    bKybdB(i) = phyvars(i).KydB*vars(i).ratio1; %Need to multiply by ratio to compare and switch sign for shoal
    bKybdB2(i) = -phyvars(i).KydB*vars(i).ratio2;
    bKyprdBpr(i) = phyvars(i).KyprdBpr*vars(i).ratio1;
    bKyprdBpr2(i) = -phyvars(i).KyprdBpr*vars(i).ratio2;
    numdays(i) = vars(i).numdays;
    Kymax(i) = vars(i).Kytop;
    bmutilBtil(i) = mean(phyvars(i).mutilBtil);
    bmutilBtil2(i) = mean(phyvars(i).mutilBtil2);
    Rk(i) = bKyprdBpr(i)./bKybdB(i);
    Bbar(i) = phyvars(i).Bbar;
    Bbar2(i) = phyvars(i).Bbar2;
    B0bar(i) = phyvars(i).B0bar;
    B0bar2(i) = phyvars(i).B0bar2;
    bioexp(i) = -((ZP)*(phyvars(i).Bbar)+(ws/H+BG/H)*B0bar(i)); %CHanged to reflect updated paper
    bioexp2(i) = -((ZP2)*(phyvars(i).Bbar2)+(ws2/H2+BG2/H2)*B0bar2(i));
    tau_bio(i) = 1./abs(phyvars(i).munetave-ZP-BG/H);
    tau_exch_Kymax(i) =(dy^2)./(vars(i).Kytop);
    tau_exch_Kymin(i) =(dy^2)./(vars(i).Kymin);
    tau_bio2(i) = 1./abs(phyvars(i).munetave2-ZP2-BG2/H2);
    tau_exch_Kymax2(i) = (dy2^2)./(vars(i).Kytop);
    tau_exch_Kymin2(i) =(dy2^2)./(vars(i).Kymin);
    total(i) = bmunetbB(i)+bmunetprBpr(i)+bKybdB(i)+bKyprdBpr(i)+bioexp(i);
    total2(i) = bmunet2bB2(i)+bmunet2prB2pr(i)+bKybdB2(i)+bKyprdBpr2(i)+bioexp2(i);
end

% load([savefile,'\','Selected PhytoBloom.mat']);
%New -24.8 and -2: 
%Bbar[-24.8] = 29.2132
%Bbar2[-24.8] = 83.3261
%Bbar[-2] = 26.8886
%Bbar2[-2days]= 80.0193