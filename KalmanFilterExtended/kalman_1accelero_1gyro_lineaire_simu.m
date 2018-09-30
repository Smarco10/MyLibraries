close all;
clear all;
clc;

format longEng

fe = 50; % fréquence d'échantillonnage
temps = 10; % temps total de mesure en seconde
t=(0:1/fe:temps-1/fe);

% Calcul des vecteur d'états à estimer (le vecteur d'état est normalement l'inconnue du système !)
vecteur_etat = zeros(temps*fe, 3); % [va ; a; b]
vecteur_etat(:, 1) = pi * 2*pi*1*sin(2*pi*1*t(:));
vecteur_etat(:, 2) = -1*pi * cos(2*pi*1*t(:)) + pi;
vecteur_etat(:, 3) = 10*t(:)+20*sin(2*pi*0.1*t(:));

% Bruit des capteurs (écart type)
bruit_capteur = zeros(2);
bruit_capteur(1) = pi*2*pi*1 * 0.03;
bruit_capteur(2) = pi*2 * 0.1;

% Calcul des paramètres mesurés
mesure = zeros(temps*fe, 2); %[va; a]
mesure(:, 1) = vecteur_etat(:, 1)+vecteur_etat(:, 3) + randn(temps*fe, 1)*bruit_capteur(1);
mesure(:, 2) = vecteur_etat(:, 2) + randn(temps*fe, 1)*bruit_capteur(2);

%[b, a] = butter(1, 5/25, 'low');
%freqz(b, a, 512, 50);
%mesure(:, 1) = filter(b, a, mesure(:, 1));
%mesure(:, 2) = filter(b, a, mesure(:, 2));

% Initialisation de Kalman
H = [1 0 1; 0 1 0];
R = [bruit_capteur(1)^2 0; 0 bruit_capteur(2)^2];
A = [1 0 0 ; 1/fe 1 0 ; 0 0 1];
Q = eye(3) * 1;
Q(1, 1) = 10;
Q(2, 2) = 0;
Q(3, 3) = 2;

X = zeros(3, 1);
X(1) = 0;%2*pi*1;
X(2) = 0;
X(3) = 0;
X_svg = zeros(temps*fe, 3);
X_svg(1, 1) = X(1);
X_svg(1, 2) = X(2);
X_svg(1, 3) = X(3);
P = zeros(3, 3);
P(1, 1) = 0;
P(2, 2) = 0;
P(3, 3) = 0;
P(1, 2) = 0;
P(2, 1) = P(1, 2);
P(2, 3) = P(1, 2);
P(3, 2) = P(1, 2);
P(1, 3) = 0;
P(3, 1) = P(1, 3);
P_svg = zeros(temps*fe, 3);
P_svg(1, 1) = P(1, 1);
P_svg(1, 2) = P(2, 2);
P_svg(1, 3) = P(3, 3);

% Calculs
for run = 2 : temps*fe
    % Prédiction
    Xp = A*X;
    Pp = A*P*A'+Q;
    
    % Mise à jour    
    K = Pp*H'*inv(R+H*Pp*H');
    P = Pp - K*H*Pp;
    X = Xp + K*(mesure(run, :)'-H*Xp);
    
    X_svg(run, 1) = X(1);
    X_svg(run, 2) = X(2);
    X_svg(run, 3) = X(3);
end

% Affichage
figure_handle=figure(1);clf;

subplot(1, 3, 1);hold on;
plot(t,mesure(:,1),'r+');
plot(t,vecteur_etat(:,1),'k');
plot(t,X_svg(:,1),'g');
legend('Mesures gyroscope','Vitesse angulaire','Estimation');
axis square;
xlabel('temps');ylabel('vitesse angulaire');
set(figure_handle,'name',' Mesures générées');

subplot(1, 3, 2);hold on;
plot(t,mesure(:,2),'r+');
plot(t,vecteur_etat(:,2),'k');
plot(t,X_svg(:,2),'g');
legend('Mesures accéléromètre', 'Angle','Estimation');
axis square;
xlabel('temps');ylabel('angle');
set(figure_handle,'name',' Mesures générées');

subplot(1, 3, 3);hold on;
plot(t,vecteur_etat(:,3),'k');
plot(t,X_svg(:,3),'g');
legend('Biais gyroscope','Estimation');
axis square;
xlabel('temps');ylabel('biais gyroscope');
set(figure_handle,'name',' Mesures générées');
