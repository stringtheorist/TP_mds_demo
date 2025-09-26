%% POROGRAMME PRINCIPAL
%% INITIALISATION =========================================================
clear;close all;clc;

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
% Tension : corde accordée sur le la-440
Note=440;       % Fréquence fondamentale [Hz]
C=2*L*Note;     % Celerité [m/s]
N0=ro*A*C^2;    % Tension [N]
Def=N0/(E*A);   % Deformation [~]
% Excitation: corde pincée d'une hauteur H en s=el 
H=L/5;          % Hauteur [m]
el=L/4;         % poistion [m]
% Domaine modal
nmax=10;        % Nombre maximal de mode considéré      
n=(1:nmax)';    % Indices modaux
Nw=nmax
kn=n*pi/L;      % Nombres d'ondes [1/m] : corde fixée aux deux extrémités
wn=C*kn;        % Pulsation [rad/s], relation de dispersion
Lamb=2*pi./kn;  % Longueur d'onde de chaque mode [m]
Per=2*pi./wn;   % Periode de chaque mode [s]
Freq=1./Per;    % Fréquence de chaque mode [Hz]
% Domaine spatial
ds=min(Lamb)/20;% Pas en espace [m]
s=(0:ds:L);     % Echantillonage spatial [m]
Ns=length(s);   % Nombre de points d'échantillonages en espace
% Domaine temporel
dt=min(Per)/20; % Pas en temps [s]
tmax=max(Per)*2;% Temps maximum de la simulation [s]
t=0:dt:tmax;    % Echantillonage temporel [s]
Nt=length(t);   % Nombre de points d'échantillonages en temps

disp(['[Nt,Ns,Nw]=[' num2str([Nt,Ns,Nw]) ']'])


Fct1(Y, s, kn);



