%% ========================================================================
%% INITIALISATION =========================================================
clear;close all;clc;
%<<FAUX>> : certaines relations utilisées ici sont fausses dans le cas 
% traité ici, nanmoins elles sont utilisée afin d'avoir des ordres de 
% grandeurs raisonnables. 

% Geometrie : section cicrculaire
L=1;            % Longueur [m]
R=0.001;        % Rayon [m]
A=pi*R^2;       % Aire [m^2]
% Materiau : acier
E=210e9;        % Module de Young [Pa]
ro=7800;        % Masse volumique [kg/m^3]
% Tension : corde accordée sur le la-440
Note=440;       % Fréquence fondamentale [Hz]
C=2*L*Note;     % Celerité [m/s] <<FAUX>>
N0=ro*A*C^2;    % Tension [N]
Def=N0/(E*A);   % Deformation [~]

% Excitation: impacte ponctuel en s=el 
V=C*2;          % Vitesse de l'impact [m/]
el=L/6;         % Position de l'impact[m]
K=15.1*pi*N0/L;  % Raideur du ressort appliqué en s=L [N/m]
% Domaine modal
nmax=10;        % Nombre maximal de mode considéré      
n=(1:nmax)';    % Indices modaux

% Domaine spectral
kmax=2*pi/R/40; % Ainsi la plus petite longueur d'onde est ~ 20xR <<FAUX>>
wmax=C*kmax;    % Plus grande longueur pulsation [rad/s] <<FAUX>>
dw=(C*pi./L)/10;% Le pas dw est ~ 1/10 de min(pulsation d'une corde fixée aux deux extrémités) <<FAUX>>
w=0:dw:wmax;    % Echantillonage en pulsation [rad/s]
Nw=length(w);   % Nombre de points d'échantillonages en w

% Domaine spatial
ds=R;         % Pas en espace [m]
s=0:ds:L;       % Echantillonage spatial [m]
Ns=length(s);   % Nombre de points d'échantillonages en espace

% Domaine temporel
dt=2*pi/(wmax/2);   % Pas en temps [s]
tmax=2*pi/(dw);% Temps maximum de la simulation [s]
t=0:dt:tmax;    % Echantillonage temporel [s]
Nt=length(t);   % Nombre de points d'échantillonages en temps
%% ========================================================================
%Rq% Dans une phase de bebeugage, il faut que [Nt,Ns,Nw] aient des valeurs 
% raisonnables (<=1000) et si possible distinctes.
disp(['[Nt,Ns,Nw]=[' num2str([Nt,Ns,Nw]) ']'])
%% ========================================================================
%% RELATION DE DISPERSION =================================================
k=sqrt(ro*A/N0)*w;
%-> visualisation de la relation de dispersion
figure(1);
plot(w,k,'LineWidth',2)
xlabel('w [rad/s]');ylabel('k [1/m]')
set(gca,'FontSize',24)

%% ========================================================================
%% EQUATION TRANSCENDENTALE ===============================================
Ka=k*L;
dKa=Ka(2)-Ka(1);
Xi=K*L/N0;
F=2*Xi*Ka./(Xi^2+Ka.^2).*cos(Ka)+sin(Ka);
%-> visualisation de la relation de dispersion
figure(2);
plot(Ka,F,'LineWidth',2)
xlabel('Ka [~]');ylabel('F')
set(gca,'FontSize',24)
% 
%% ========================================================================
%% RECHERCHE DES RACINES DE L'EQUATION TRANSCENDENTALE ====================
% On cherche d'abord numériquement une bonne approximation de ces racines
signF=sign(F);
diffsignF=diff(signF);
[iKa]=find(abs(diffsignF)==2)+1;
nmax=length(iKa);
%-> visualisation de la relation de dispersion
figure(3);hold on
plot(Ka,F,'LineWidth',2)
plot(Ka,signF,'LineWidth',2)
plot(Ka(1:Nw-1),diffsignF,'LineWidth',2)
plot(Ka(iKa),F(iKa),'o','MarkerSize',8,'LineWidth',2)
xlabel('Ka [~]');ylabel('F')
set(gca,'FontSize',24)

% On utilise ensuite un algorithme de cherche de zéros de fonctions.  avec 
% Le domaine de recherche est situé dans un voisinage de l'approximation 
% trouvée précédement. 
%-> La fonction quadratique à minimiser
fctFquad=@(Ka) (2*Xi*Ka./(Xi^2+Ka.^2).*cos(Ka)+sin(Ka)).^2;
%-> La minimisation
for n=1:nmax
    Kan(n)=fminbnd(@(KK) fctFquad(KK),Ka(iKa(n))-2*dKa,Ka(iKa(n))+2*dKa);    
end
%-> visualisation de la relation de dispersion
figure(4);hold on
plot(Ka,F,'LineWidth',2)
plot(Kan,2*Xi*Kan./(Xi^2+Kan.^2).*cos(Kan)+sin(Kan),'o','MarkerSize',8,'LineWidth',2)
xlabel('Ka [~]');ylabel('F')
set(gca,'FontSize',24)
%% ========================================================================

%% ========================================================================
%% FREQUENCE PROPRES ET RE-DIMENSIONNALISATION DES DOMAINES ===============
n=(1:nmax)';    % Indices modaux
kn=Kan'/L;      % Nombres d'ondes [1/m]
wn=C*kn;        % Pulsation [rad/s], relation de dispersion
Lamb=2*pi./kn;  % Longueur d'onde de chaque mode [m]
Per=2*pi./wn;   % Periode de chaque mode [s]
Freq=1./Per;    % Fréquence de chaque mode [Hz]
% Domaine spatial
ds=min(Lamb)/20;% Pas en espace [m]
s=(0:ds:L);     % Echantillonage spatial [m]
% Domaine temporel
dt=min(Per)/2;  % Pas en temps [s]
tmax=max(Per)*2;% Temps maximum de la simulation [s]
t=0:dt:tmax;    % Echantillonage temporel [s]
%% ========================================================================

%% ========================================================================
%% ANALYSE MODALE =========================================================
% Modes propres
for in=1:nmax
    % Y_ij, avec i=>n et j=>s
Y(in,:)=cos(kn(in).*s)+K/N0*sin(kn(in).*s)./kn(in);    % Y_ij, avec i=>n et j=>s
% Norme au carré de chaque mode (obtenu via calcul formel) 
normY2(in,:)=(2*(K^2.*kn(in)*L+K*N0+kn(in)*L*N0^2)-2*K*N0*cos(2*kn(in)*L)+(N0^2-K^2).*sin(2*kn(in)*L))./(4*kn(in)*N0^2);
end
%-> visualisation de quelques modes propres
figure(5);
plot(s,Y([1:3 nmax],:),'LineWidth',2)
xlabel('s [m]')
legend('n=1','n=2','n=3','n=nmax')
set(gca,'FontSize',24)

%% Amplitude modale
an=zeros(size(n));
bn=L*V*(cos(kn.*el)+K/N0*sin(kn.*el)./kn)./normY2;
%-> visualisation des amplitudes modales an
figure(6);
stem(Freq,abs(bn),'LineWidth',2)
xlabel('fn [Hz]')
ylabel('|bn| [m]')
set(gca,'FontSize',24)

%% Fonction en temps
T=an.*cos(wn.*t)+bn.*sin(wn.*t)./wn;    % T_ij, avec i=>n et j=>t
%-> visualisation de T(t) pour quelques modes
figure(7);
plot(t,T([1:3 nmax],:),'LineWidth',2)
xlabel('t [s]')
legend('n=1','n=2','n=3','n=nmax')
set(gca,'FontSize',24)

%% Solution générale
u=Y'*T;    % u_ij, avec i=>s et j=>t
%-> visualisation de u(s,t) à divers instants
figure(8);
plot(s,u(:,[1 10 20]),'LineWidth',2);
xlabel('s [m]');ylabel('u(s,t) [m]');
legend(['t=' num2str(t(1)) ],['t=' num2str(t(10)) ],['t=' num2str(t(20)) ])
axis equal
set(gca,'FontSize',24)
%-> visualisation de u(s,t) en divers point de la corde
figure(8);
plot(t,u([1 10 20],:),'LineWidth',2);
xlabel('t [m]');ylabel('u(s,t) [m]');
legend(['s=' num2str(s(1)) ],['s=' num2str(s(10)) ],['s=' num2str(s(20)) ])
set(gca,'FontSize',24)
%-> visualisation de u(s,t) au cours du temps
figure(9);
for j=1:length(t)
    plot(s,u(:,j),'LineWidth',2);
    xlabel('s [m]');ylabel('u(s,t) [m]');
    axis equal;axis([0,L,-1*1*max(max(abs(u))),1*1*max(max(abs(u)))])
    set(gca,'FontSize',24);
    pause(0.1)
end