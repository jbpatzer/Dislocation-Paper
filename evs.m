function [mu,mu1,chi,psi,dpsi,dchi] = evs(wl,D,L,plotflag)
%% Find the eigenvalues and eigenfunctions of PSP and PSP'
% wl = wavelength
% D = diffusion parameter
% L = periodic solution
% plotflag = true for plotting

N = length(L);
x = linspace(0,wl,N+1); x = x(2:end); dx = x(2)-x(1);
u = exp(L);

% FFT
k = -N/2:N/2-1;
FTphi = ifftshift(sinc(k/wl).*(-1).^k);
J0 = toeplitz(ifftshift((ifft(fft([1 zeros(1,N-1)]).*FTphi))));

% Spectral Differentiation matrices
[~, Mz] = fourdif(N,1);
[~, M] = fourdif(N,2);
M = 4*pi^2*M/wl^2; Mz = 2*pi*Mz/wl;
Nu = ifftshift((ifft(fft(u).*FTphi)));

%% solve PSP
S = D*M+spdiags(1-Nu',0,N,N)-u'.*J0';
[psi,mu] = eig(S); mu = -diag(mu)/D;
[~,ind] = sort(real(mu),'ascend');
mu = mu(ind);
psi = psi(:,ind);
psi = psi.*sign(real(psi(N/2+10,:)));

%% solve PSP'
S1 = D*M+spdiags(1-Nu',0,N,N)-J0'.*u;
[chi,mu1] = eig(S1); mu1 = -diag(mu1)/D;
[~,ind] = sort(real(mu1),'ascend');
mu1 = mu1(ind);
chi = chi(:,ind);
chi = chi.*sign(real(chi(N/2+10,:)));

%% Make sure eigenvalues are ordered correctly and normalize
flag = true;
for k = 1:N
    if flag&&(k<N)&&(abs(imag(mu(k)))>1e-12)&&(abs(imag(mu(k+1)))>1e-12)&&(abs(mu(k)-mu1(k+1))>1e-12)
        mu1(k:k+1) = conj(mu1(k:k+1)); chi(:,k:k+1) = conj(chi(:,k:k+1));
        flag = false;
    else
        flag = true;
    end
    % Normalize eigenfunctions
    cc = dx*chi(:,k)'*chi(:,k);
    pp = dx*psi(:,k)'*psi(:,k);
    cp = dx*chi(:,k)'*psi(:,k);

    fac1 = sqrt(sqrt(pp/cc)); acp = abs(cp);
    b = cp*fac1/acp^(3/2); a = 1/fac1/sqrt(acp);

    psi(:,k) = a*psi(:,k); chi(:,k) = b*chi(:,k);

end
dpsi = Mz*psi; dchi = Mz*chi;

if plotflag
    figure(2),clf
    semilogy(1:length(mu)-1,abs(mu(2:end)),'.k')
    grid on
    ylabel('|\mu|')
    grid on
    title(sprintf('D = %1.2e, \\lambda = %1.3f',D,wl))
    set(gca,'FontSize',18)

    figure(3),clf
    subplot(3,1,1)
    plot(x,psi(:,1:5),'LineWidth',1.2)
    grid on
    xlabel('x')
    ylabel('\psi')
    title(sprintf('D = %1.2e, \\lambda = %1.3f, \\mu = %1.2e, %1.2e, %1.2e, %1.2e, %1.2e...',D,wl,mu(1),mu(2),mu(3),mu(4),mu(5)))
    legend('0','1','2','3','4')
    set(gca,'FontSize',18)
    subplot(3,1,3)
    plot(x,chi(:,1:5),'LineWidth',1.2)
    grid on
    xlabel('x')
    ylabel('\chi')
    set(gca,'FontSize',18)

    subplot(3,1,2)
    plot(x,dpsi(:,1:5),'LineWidth',1.2)
    grid on
    xlabel('x')
    ylabel('\psi''')
    set(gca,'FontSize',18)

    figure(4),clf
    for k = 1:9
        subplot(9,1,k)
        plot(x,psi(:,k),'k','linewidth',1.2)
        grid on
        ylabel(['\psi_' num2str(k-1)])
        set(gca,'FontSize',18)
    end
    xlabel('x')
    subplot(9,1,1)
    title(sprintf('D = %1.2e, \\lambda = %1.3f',D,wl))
    set(gca,'FontSize',18)

    figure(5),clf
    for k = 1:9
        subplot(9,1,k)
        plot(x,chi(:,k),'k','linewidth',1.2)
        grid on
        ylabel(['\chi_' num2str(k-1)])
        set(gca,'FontSize',18)
    end
    xlabel('x')
    subplot(9,1,1)
    title(sprintf('D = %1.2e, \\lambda = %1.3f',D,wl))

    figure(6),clf
    for k = 1:9
        subplot(9,1,k)
        plot(x,dpsi(:,k),'k','linewidth',1.2)
        grid on
        ylabel(['\psi''_' num2str(k-1)])
        set(gca,'FontSize',18)
    end
    xlabel('x')
    subplot(9,1,1)
    title(sprintf('D = %1.2e, \\lambda = %1.3f',D,wl))
    set(gca,'FontSize',18)
end