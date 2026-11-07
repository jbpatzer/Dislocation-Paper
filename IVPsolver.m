function IVPsolver(N,wl,D,tend,dtout,moviefile,movietime,umax,plotflag,contflag,showdiff,Ninitial,delta)
%% Solve the initial value problem (IVP), formulated in terms of L = log(u)

% N = number of grid points
% wl  = length of computational domain (solution has this period)
% D = diffusion parameter
% tend - Solve for 0 < t < tend
% moviefile - name of movie file - empty for no movie
% movietime = length of movie is seconds
% umax - maximum value of y-axis in output plot
% plotflag - plot if true
% contflag - continue from crashed run if true
% showdiff - also plot difference from unperturbed solution if true
% Ninitial - number of initial peaks
% delta - dislocation length

format short e
format compact
t0 = toc;
Db = Dbar(2^9,wl/Ninitial,D,[]);

if contflag
    load('cont','t','tout','x','L','uold','xpeaks','k1');
    xold = x; u = exp(L); 
else
    xold = []; uold = []; 
    % Initial condition

    wl0 = wl/Ninitial; % underlying wavelength
    % Find F_p and repeat periodically to find unperturbed solution
    [L,x0] = periodic(2^10,wl0,D,[],false);
    pp = csape(x0, L, 'periodic');
    
    %Dislocate
    wl = wl+delta;
    x = linspace(-0.5*wl,0.5*wl,N+1);
    L = displace(x,pp,wl0,delta);

    u = exp(L);
    u = u(2:end); L = L(2:end);
    x = x(2:end);

    t = 0; tout = 0; tout2=0;
    u0 = exp(ppval(pp,mod(x,wl0))); % Unperturbed solution
    du0 = u-u0;

end

dx = x(2)-x(1);

% precompute FT of kernel
k = -N/2:N/2-1;
FTphi = ifftshift(sinc(k/wl).*(-1).^k);

if ~isempty(moviefile)
    set(gca,'nextplot','replacechildren');
    v = VideoWriter(moviefile,'MPEG-4');
    Nout = ceil(tend/dtout);
    v.FrameRate = max(5,min(48,ceil(Nout/movietime)));
    open(v);
else
    v = [];
end

if plotflag
    output(x,L,wl,D,t,moviefile,v,umax,xold,uold,u0,du0,delta,wl0,Db,showdiff)
end

%% Set up fourth order accurate differentiation matrices
M = spdiags([-ones(N,1) 16*ones(N,1) -30*ones(N,1) 16*ones(N,1) -ones(N,1)],-2:2,N,N);
M(1,N-1) = -1; M(1,N) = 16; 
M(2,N) = -1;
M(N,1) = 16; M(N,2) = -1; 
M(N-1,1) = -1;
Mz = spdiags([ones(N,1) -8*ones(N,1) 8*ones(N,1) -ones(N,1)],[-2 -1 1 2],N,N);
Mz(1,N-1) = 1; Mz(1,N) = -8; 
Mz(2,N) = 1;
Mz(N,1) = 8; Mz(N,2) = -1; 
Mz(N-1,1) = -1;
M = M/dx^2/12;
Mz = Mz/dx/12;

difftol = 0.1;
dtmax = min(0.1,dx^2/D/8);
dt = min(0.01,dtmax); dtmin = 1e-4; 

while t<tend
    %Explicit midpoint method
    q = RHS(L,FTphi)+D*(M*L'+(Mz*L').^2);
    Lhalf =  L+0.5*dt*q';
    q = RHS(Lhalf,FTphi)+D*(M*Lhalf'+(Mz*Lhalf').^2);
    Lnew = L+dt*q';
    unew = exp(Lnew);
    dumax = max(abs(unew-u)); 

    if (dumax>=difftol)&&(dt>dtmin)
        dt = max(dtmin,0.9*difftol*dt/dumax);
    else
        t = t+dt; tout = tout+dt; tout2 = tout2+dt; L = Lnew; u = unew;

        if plotflag&&(tout>dtout)
            tout = tout-dtout;
            output(x,L,wl,D,t,moviefile,v,umax,xold,uold,u0,du0,delta,wl0,Db,showdiff)
            save cont x t L tout uold u0
        end

        if (dumax<0.8*difftol)&&(t>10)
            dt = min(dtmax,1.1*dt);
        end
    end
end

if plotflag
    output(x,L,wl,D,t,moviefile,v,umax,xold,uold,u0,du0,delta,wl0,Db,showdiff)
end

if ~isempty(moviefile)
    close(v);
end

fprintf('Elapsed time = %3.2f minutes.\n',(toc-t0)/60)

%%
function RHS = RHS(L,FTphi)
RHS = (1-ifftshift(real(ifft(fft(exp(L)) .* FTphi))))'; 

%%
function L = displace(x,pp,wl0,delta)
L = zeros(size(x));
for k = 1:length(x)
    if abs(x(k))<0.5*delta
        L(k) = ppval(pp,0);
    else
        L(k) = ppval(pp,mod(x(k)-sign(x(k))*0.5*delta,wl0));
    end
end

%%
function output(x,L,wl,D,t,moviefile,v,umax,xold,uold,u0,du0,delta,wl0,Db,showdiff)
persistent dumax
figure(1)
u = exp(L);
if isempty(umax)
    umax = 1.1*max(u);
else
    umax = max(umax,1.1*max(u));
end
figure(1),clf
if showdiff
    subplot(2,1,1)
    du = u-u0;
    plot(x,du0,'k--',x,du,'k','LineWidth',1.2);
    xlim([-0.5*wl 0.5*wl])
    if delta>0
        ylim([min(-umax*10*delta,1.1*min(du)) max(umax*10*delta,1.1*max(du))])
    end
    ylabel('$u(x,t;\delta)-u(x,0;0)$','interpreter','latex')
    grid on
    set(gca,'FontSize',18)
    subplot(2,1,2)
end

plot(x,u,'k','LineWidth',1.2);
axis([-0.5*wl 0.5*wl 0 umax])
title(sprintf('$t = %6.0f, D = 10^{%.g}, \\lambda = %1.4f, \\delta = 10^{%.g}, \\bar{D} = %1.5f$',t,log10(D),wl0,log10(delta),Db),'interpreter','latex')
xlabel('x','interpreter','latex')
ylabel('$u(x,t)$','interpreter','latex')
grid on
set(gca,'FontSize',18)
drawnow
if ~isempty(moviefile)
    writeVideo(v,getframe(gcf));
end
if ~isempty(uold)
    figure(4),clf
    du = uold-spline(x,u,xold);
    plot(xold,du,'k','linewidth',1.2)
    grid on
    xlabel('x')
    ylabel('du')
    if isempty(dumax)
        dumax = max(du);
    end
    ylim([-dumax dumax])
    drawnow
else
    dumax = [];
end