%Biology and Exchange Timescales
close all

depamu = 1/H*sum(dz.*munetm);
depamu2 = 1/H2*sum(dz.*munetm2);

tau_bio = 1./abs(depamu-ZP-BG/H);
tau_exch =(dy^2)./(Kym(1,:));
tau_bio2 = 1./abs(depamu2-ZP2-BG2/H2);
tau_exch2 = (dy2^2)./(Kym(1,:));

figure()
p1 = plot(times/24/3600,tau_bio/24/3600,'LineWidth',3);
ylabel('Timescale (d)')
xlabel('Time (d)')
hold on
p2 = plot(times/24/3600,tau_exch/24/3600);
%p3 = plot(times/24/3600,tau_bio./tau_exch);
hold off
leg1 = legend('$\tau_{bio}=\frac{1}{abs(\langle \mu_{net,1}\rangle -ZP_1-\frac{BG_1}{H_1})}$','$\tau_{exch}=\frac{dy_1^2}{Ky}$'); %,'$\frac{\tau_{bio}}{\tau_{exch}}$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',12);
title(['Channel Biology vs Exchange Timescales Kymax = ', num2str(Ky_range(phyvars(run_num).Kytest)), ', Period = ', num2str(num_range(phyvars(run_num).daystest)*2)])

set(gca, 'SortMethod', 'depth')
p1.ZData = ones(size(p1.XData));
p2.ZData = zeros(size(p2.XData));

figure()
p1 = plot(times/24/3600,tau_bio2/24/3600,'LineWidth',3);
ylabel('Timescale (d)')
xlabel('Time (d)')
hold on
p2 = plot(times/24/3600,tau_exch2/24/3600);
%p3 = plot(times/24/3600,tau_bio2./tau_exch2);
hold off
leg1 = legend('$\tau_{bio}=\frac{1}{abs(\langle \mu_{net,2}\rangle -ZP_2-\frac{BG_2}{H_2})}$','$\tau_{exch}=\frac{dy_2^2}{Ky}$'); %,'$\frac{\tau_{bio}}{\tau_{exch}}$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',12);
title(['Shoal Biology vs Exchange Timescales Kymax = ', num2str(Ky_range(phyvars(run_num).Kytest)), ', Period = ', num2str(num_range(phyvars(run_num).daystest)*2)])


set(gca, 'SortMethod', 'depth')
p1.ZData = ones(size(p1.XData));
p2.ZData = zeros(size(p2.XData));


%WAY CALCULATED IN PAPER
tbio_paper = mean(tau_bio(teqstart:tend))/86400;
tbio_paper2 = mean(tau_bio2(teqstart:tend))/86400;

%% Another option exchange scale

tau_bio = 1./abs(depamu-ZP-BG/H);
tau_exch = depaBio.*(dy^2)./(Kym(1,:).*(depaBio2-depaBio));
tau_bio2 = 1./abs(depamu2-ZP2-BG2/H2);
tau_exch2 = depaBio2.*(dy2^2)./(Kym(1,:).*(depaBio2-depaBio));

figure()
p1 = plot(times/24/3600,tau_bio/24/3600,'LineWidth',3);
ylabel('Timescale (d)')
xlabel('Time (d)')
hold on
p2 = plot(times/24/3600,tau_exch/24/3600);
%p3 = plot(times/24/3600,tau_bio./tau_exch);
hold off
leg1 = legend('$\tau_{bio}=\frac{1}{abs(\langle \mu_{net,1}\rangle -ZP_1-\frac{BG_1}{H_1})}$','$\tau_{exch}=\frac{\langle B_{1}\rangle dy_1^2}{Ky(\langle B_2\rangle-\langle B_1\rangle)}$'); %,'$\frac{\tau_{bio}}{\tau_{exch}}$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',12);
title(['Channel Biology vs Exchange Timescales Kymax = ', num2str(Ky_range(phyvars(run_num).Kytest)), ', Period = ', num2str(num_range(phyvars(run_num).daystest)*2)])

set(gca, 'SortMethod', 'depth')
p1.ZData = ones(size(p1.XData));
p2.ZData = zeros(size(p2.XData));

figure()
p1 = plot(times/24/3600,tau_bio2/24/3600,'LineWidth',3);
ylabel('Timescale (d)')
xlabel('Time (d)')
hold on
p2 = plot(times/24/3600,tau_exch2/24/3600);
%p3 = plot(times/24/3600,tau_bio2./tau_exch2);
hold off
leg1 = legend('$\tau_{bio}=\frac{1}{abs(\langle \mu_{net,2}\rangle -ZP_2-\frac{BG_2}{H_2})}$','$\tau_{exch}=\frac{\langle B_{2}\rangle dy_2^2}{Ky(\langle B_2\rangle-\langle B_1\rangle)}$'); %,'$\frac{\tau_{bio}}{\tau_{exch}}$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',12);
title(['Shoal Biology vs Exchange Timescales Kymax = ', num2str(Ky_range(phyvars(run_num).Kytest)), ', Period = ', num2str(num_range(phyvars(run_num).daystest)*2)])


set(gca, 'SortMethod', 'depth')
p1.ZData = ones(size(p1.XData));
p2.ZData = zeros(size(p2.XData));