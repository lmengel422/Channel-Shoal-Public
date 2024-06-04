% ********************************************************************
%       Pre-allocates all matrices for speed for watercolumn model
%       Fall 2019 Lily Engel
% ********************************************************************
if isave == 1
    mm = floor(M/isave)-1; %Size of what will be saved in isave. Need to be minus 1 if isave 1 since loop start from 2
else
    mm = floor(M/isave);
end

% *******************************************************************
%  Pre-define matrices used in watercolumn.m
% *******************************************************************
t = zeros(1,M);
times = zeros(1,mm);

if first_water_column == 1
    Um = zeros(N,mm);
    Cm = zeros(N,mm);
    Biom = zeros(N,mm);
    Q2m = zeros(N,mm);
    Q2Lm = zeros(N,mm);
    rhom = zeros(N,mm);
    Lm = zeros(N,mm);
    nu_tm = zeros(N,mm);
    Kzm = zeros(N,mm);
    Kqm = zeros(N,mm);
    N_BVm = zeros(N,mm);
    Cplot = zeros(N,mm);
end

if second_water_column == 1
    Um2 = zeros(N2,mm);
    Cm2 = zeros(N2,mm);
    Biom2 = zeros(N2,mm);
    Q2m2 = zeros(N2,mm);
    Q2Lm2 = zeros(N2,mm);
    rhom2 = zeros(N2,mm);
    Lm2 = zeros(N2,mm);
    nu_tm2 = zeros(N2,mm);
    Kzm2 = zeros(N2,mm);
    Kqm2 = zeros(N2,mm);
    N_BVm2 = zeros(N2,mm);
    Cplot2 = zeros(N2,mm);
end

if first_water_column == 1 && second_water_column == 1
    Kym = zeros(N2,mm);
    KySm = zeros(N2,mm);
    diffUm = zeros(N2,mm);
    diffCm = zeros(N2,mm);
    diffBiom = zeros(N2,mm);
    diffQ2m = zeros(N2,mm);
    diffQ2Lm = zeros(N2,mm);
    diffrhom = zeros(N2,mm);
    diffLm = zeros(N2,mm);
    diffnu_tm = zeros(N2,mm);
    diffKzm = zeros(N2,mm);
    diffKqm = zeros(N2,mm);
    diffN_BVm = zeros(N2,mm);
end

% *******************************************************************
%  Pre-define matrices used in wc_advance.m
% *******************************************************************
if first_water_column == 1
    Px = zeros(1,N);
    Sm = zeros(1,N);
    Sh = zeros(1,N);
    rho = zeros(1,N);
    N_BV = zeros(1,N);
    Q = zeros(1,N);
    L = zeros(1,N);
    Kq = zeros(1,N);
    nu_t = zeros(1,N);
    Kz = zeros(1,N);
end

if second_water_column == 1
    Px2 = zeros(1,N2);
    Sm2 = zeros(1,N2);
    Sh2 = zeros(1,N2);
    rho2 = zeros(1,N2);
    N_BV2 = zeros(1,N2);
    Q_2 = zeros(1,N2);
    L2 = zeros(1,N2);
    Kq2 = zeros(1,N2);
    nu_t2 = zeros(1,N2);
    Kz2 = zeros(1,N2);
end

% *******************************************************************
%  Pre-define matrices used in wc_setup.m
% *******************************************************************
if first_water_column == 1
    aU = zeros(1,N);
    bU = zeros(1,N);
    cU = zeros(1,N);
    dU = zeros(1,N);
    aC = zeros(1,N);
    bC = zeros(1,N);
    cC = zeros(1,N);
    dC = zeros(1,N);
    aBio = zeros(1,N);
    bBio = zeros(1,N);
    cBio = zeros(1,N);
    dBio = zeros(1,N);
    aQ2 = zeros(1,N);
    bQ2 = zeros(1,N);
    cQ2 = zeros(1,N);
    dQ2 = zeros(1,N);
    aQ2L = zeros(1,N);
    bQ2L = zeros(1,N);
    cQ2L = zeros(1,N);
    dQ2L = zeros(1,N);
end

if second_water_column == 1
    aU2 = zeros(1,N2);
    bU2 = zeros(1,N2);
    cU2 = zeros(1,N2);
    dU2 = zeros(1,N2);
    aC2 = zeros(1,N2);
    bC2 = zeros(1,N2);
    cC2 = zeros(1,N2);
    dC2 = zeros(1,N2);
    aBio2 = zeros(1,N2);
    bBio2 = zeros(1,N2);
    cBio2 = zeros(1,N2);
    dBio2 = zeros(1,N2);
    aQ22 = zeros(1,N2);
    bQ22 = zeros(1,N2);
    cQ22 = zeros(1,N2);
    dQ22 = zeros(1,N2);
    aQ2L2 = zeros(1,N2);
    bQ2L2 = zeros(1,N2);
    cQ2L2 = zeros(1,N2);
    dQ2L2 = zeros(1,N2);
end

% *******************************************************************
%  Set initial conditions for Ky and KyS- outside flag so still runs
% *******************************************************************
Ky = zeros(1,N2);
KyS = zeros(1,N2);
exchange = zeros(2,N2);
