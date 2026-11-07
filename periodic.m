function [L,x] = periodic(N,wl,D,L0,plotflag)
%% Find the solution of the steady periodic problem
%% Solves for L = log(u)

% N = number of grid points
% wl = wavelength 
% D = diffusion parameter
% L0 = initial guess of solution - will be interpolated onto any new grid - use [] for built in guess
% plotflag = true for plotting

x = linspace(0,wl,N+1); x = x(2:end); 

% FFT
k = -N/2:N/2-1;
FTphi = ifftshift(sinc(k/wl).*(-1).^k);
J0 = toeplitz(ifftshift((ifft(fft([1 zeros(1,N-1)]).*FTphi))));

% initial guess
D0 = -wl^3*sin(pi/wl)/4/pi^3;
k1 = -wl*sin(pi/wl)/pi; k2 = -wl*sin(2*pi/wl)/2/pi;
g = (k1+k2)/(k2-4*k1)+2-2*k1;
A = 2*sqrt((D0-D)/D0/g);
if A<0.5
    u = 1+A*cos(2*pi*(x-wl/2)/wl);
    L0 = log(u);
elseif isempty(L0)
    L0 = log(3)-100*(x-0.5*wl).^2/wl^2;
elseif length(L0)~=length(x)
    L0 = spline(linspace(0,wl,length(L0)),L0,x*wl/x(end));
end

% Spectral differentiation matrices
[~, Mz] = fourdif(N,1);
[~, M] = fourdif(N,2);
M = 4*pi^2*M/wl^2; Mz = 2*pi*Mz/wl;

options = optimoptions('fsolve','MaxFunEval',250,'FunctionTolerance',1e-14,...
    'SpecifyObjectiveGradient',true,'display','off');

L = fsolve(@(L) eqn(L,FTphi,M,Mz,D,J0),L0,options);
u = exp(L); 

if  plotflag
    figure(1), clf
    subplot(1,2,1)
    plot(x,L,'k', 'LineWidth',1.2);
    xlim([0 wl])
    title(sprintf('wavelength = %1.3f, D = %1.2e',wl,D))
    xlabel('x'), ylabel('L')
    grid on
    set(gca,'FontSize',18)
    subplot(1,2,2)
    plot(x,u,'k', 'LineWidth',1.2);
    xlim([0 wl])
    title(sprintf('wavelength = %1.4f, D = %1.2e',wl,D))
    xlabel('x'), ylabel('u')
    ylim([0 1.1*max(u)])
    grid on
    set(gca,'FontSize',18)
    drawnow
end

%%
function [f, J] = eqn(L,FTphi,M,Mz,D,J0)
% Nonlinear algebraic equations from discretization
u = exp(L);
Nu = ifftshift((ifft(fft(u).*FTphi)));
Lz = Mz*L';
f = D*(M*L'+Lz.^2)+1-Nu';

if nargout>1
    J = D*(M+2*Mz.*Lz)-J0.*u;
end