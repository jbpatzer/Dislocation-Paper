function Db = Dbar(N,wl,D,L0)
%% Calculate Dbar
% N = number of grid points
% wl = wavelength 
% D = diffusion parameter
% L0 = initial guess of solution - will be interpolated onto any new grid - use [] for built in guess

%% Find periodic solution
[L,x] = periodic(N,wl,D,L0,false);

%% Find eigenvalues and eigensolutions
[mu,~,chi,psi,dpsi,dchi] = evs(wl,D,L,false);
Fp = exp(L);
dx = x(2)-x(1);

%% alpha_r
ar = dx*chi'*dpsi(:,1);

%% beta_r
k = -N/2:N/2-1;
FTphi2 = -wl*ifftshift((cos(pi*k/wl)-sinc(k/wl)).*(-1).^k./(2*pi*1i*k));
FTphi2(1) = 0;
Npsi = ifftshift(ifft(fft(psi(:,1).').* FTphi2));

br = dx*chi'*(Npsi.*Fp).';

%% Gp
fac = (br/D+2*ar)./mu;
Gp = psi(:,2:end)*fac(2:end);

%% alphabar
NGp = ifftshift(ifft(fft(Gp.').* FTphi2));
a = dx*sum(chi(:,1).'.*NGp.*Fp);

%% betabar
FTphi3 = ifftshift((0.25*sinc(k/wl) + wl^2*(cos(pi*k/wl)-sinc(k/wl))/2/pi^2./k.^2).*(-1).^k);
FTphi3(1) = 1/12;

Npsi = ifftshift(ifft(fft(psi(:,1).').* FTphi3));
b = 0.5*dx*sum(chi(:,1).'.*Npsi.*Fp);

%% gammabar
g = -2*dx*sum(Gp.*dchi(:,1));

Db = D*(1+g)+a-b;