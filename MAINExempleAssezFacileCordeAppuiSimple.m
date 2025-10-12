%% ========================================================================
%% INITIALISATION =========================================================
clear;close all;clc;
% Geometrie : section cicrculaire
L=1;            % Longueur [m]
R=0.001;        % Rayon [m]
A=pi*R^2;       % Aire [m^2]
% Materiau : acier
E=210e9;        % Module de Young [Pa]
ro=7800;        % Masse volumique [kg/m^3]
% Tension : corde accordee sur le la-440
Note=440;       % Frequence fondamentale [Hz]
C=2*L*Note;     % Celerite [m/s]
N0=ro*A*C^2;    % Tension [N]
Def=N0/(E*A);   % Deformation [~]
% Excitation: corde pincee d'une hauteur H en s=el 
H=L/5;          % Hauteur [m]
el=L/4;         % poistion [m]
% Domaine modal
nmax=10;        % Nombre maximal de mode considere      
n=(1:nmax)';    % Indices modaux
Nw=nmax;
kn=n*pi/L;      % Nombres d'ondes [1/m] : corde fixee aux deux extremites
wn=C*kn;        % Pulsation [rad/s], relation de dispersion
Lamb=2*pi./kn;  % Longueur d'onde de chaque mode [m]
Per=2*pi./wn;   % Periode de chaque mode [s]
Freq=1./Per;    % Frequence de chaque mode [Hz]
% Domaine spatial
ds=min(Lamb)/20;% Pas en espace [m]
s=(0:ds:L);     % Echantillonage spatial [m]
Ns=length(s);   % Nombre de points d'echantillonages en espace
% Domaine temporel
dt=min(Per)/20; % Pas en temps [s]
tmax=max(Per)*2;% Temps maximum de la simulation [s]
t=0:dt:tmax;    % Echantillonage temporel [s]
Nt=length(t);   % Nombre de points d'echantillonages en temps
% %% ========================================================================
%Rq% Dans une phase de bebeugage, il faut que [Nt,Ns,Nw] aient des valeurs 
% raisonnables (<=1000) et si possible distinctes.
disp(['[Nt,Ns,Nw]=[' num2str([Nt,Ns,Nw]) ']'])

%% ========================================================================
%% ANALYSE MODALE =========================================================
% Modes propres
for in=1:Nw
    % Y_ij, avec i=>n et j=>s
    Y(in,:)=sin(kn(in)*s); 
end
%=> visualisation de quelques modes propres
figure(1);
plot(s,Y([1:3 nmax],:),'LineWidth',2)
xlabel('s [m]')
legend('n=1','n=2','n=3','n=nmax')
set(gca,'FontSize',24)

% Amplitude modale
an=2*H./(n*pi)*L/(L-el).*sin(kn*el)./(kn*el);
bn=zeros(size(n));
%=> visualisation des amplitudes modales an
figure(2);
stem(wn,abs(an),'LineWidth',2)
xlabel('wn [rad/s]')
ylabel('|an| [m]')
set(gca,'FontSize',24)

% Fonction en temps
for in=1:nmax
    % T_ij, avec i=>n et j=>t
    T(in,:)=an(in)*cos(wn(in)*t)+bn(in)*sin(wn(in)*t)/wn(in); 
end
%=> visualisation de T(t) pour quelques modes
figure(3);
plot(t,T([1:3 nmax],:),'LineWidth',2)
xlabel('t [s]')
legend('n=1','n=2','n=3','n=nmax')
set(gca,'FontSize',24)
%% ========================================================================

%% ========================================================================
%% SOLUTION GENERALE ======================================================
u=Y'*T;    % u_ij, avec i=>s et j=>t
%-> visualisation de u(s,t) a divers instants
figure(4);subplot(1,2,1)
plot(s,u(:,[1 10 20]),'LineWidth',2);
xlabel('s [m]');ylabel('u(s,t) [m]');
legend(['t=' num2str(t(1)) ],['t=' num2str(t(10)) ],['t=' num2str(t(20)) ])
axis equal
set(gca,'FontSize',24)
%-> visualisation de u(s,t) en divers point de la corde
figure(4);subplot(1,2,2)
plot(t,u([1 10 20],:),'LineWidth',2);
xlabel('t [m]');ylabel('u(s,t) [m]');
legend(['s=' num2str(s(1)) ],['s=' num2str(s(10)) ],['s=' num2str(s(20)) ])
set(gca,'FontSize',24)
%-> visualisation de u(s,t) au cours du temps
figure(5);
for j=1:Nt
    plot(s,u(:,j),'k','LineWidth',2);hold on
    plot(s([1 10 20]),u([1 10 20],j),'o','MarkerSize',8,'LineWidth',2)
    hold off
    xlabel('s [m]');ylabel('u(s,t) [m]');
    axis equal;axis([0,L,-H,H])
    set(gca,'FontSize',24);
    pause(0.1)
end
%% ========================================================================
