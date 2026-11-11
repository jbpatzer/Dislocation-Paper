D = 1e-4; N = 2^9; M = 150;
wl = [0.52 0.53 0.54 0.55];
figure(1),clf
for k = 1:4
    subplot(3,2,k)
    gamma = linspace(0,pi/wl(k),M);
    evs = Bloch(D, wl(k), N, gamma, false);
    Db = real(Dbar(N,wl(k),D,[]));
    plot(gamma,real(evs(1,:)),'k',-gamma,real(evs(1,:)),'k-',...
        gamma, Db*gamma.^2,'k--',-gamma,Db*gamma.^2,'k--','linewidth',1.4)
    if k>2
        xlabel('$\gamma$','interpreter','latex')
    end
    if (k==1)||(k==3)
        ylabel('$\sigma$','interpreter','latex')
    end
    title(sprintf('$\\lambda$ = %.3g', wl(k)),'interpreter','latex')
    set(gca,'FontSize',18)
    grid on
    ylim([-pi/wl(k) pi/wl(k) -0.2 0.2])
    drawnow
end

wl = linspace(wl(1),wl(end),M);
Db = zeros(1,M);
for k = 1:M
    Db(k) = real(Dbar(N,wl(k),D,[]));
end

subplot(3,1,3)
plot(wl,Db,'k','linewidth',1.4)
grid on
xlabel('$\lambda$','interpreter','latex')
xlabel('$\bar{D}$','interpreter','latex')
set(gca,'FontSize',18)
xlim([wl(1) wl(end)])