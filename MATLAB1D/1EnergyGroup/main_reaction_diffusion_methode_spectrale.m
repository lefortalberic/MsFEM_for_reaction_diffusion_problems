% =====================================================
%
% principal_reaction_diffusion;
%
% une routine qui calcul pour un pas H donné l'erreur entre la solution approchée par la méthode de résolution de notre choix
% et la solution exacte, calculée avec les éléments P1 classiques de pas de
% maillage très petit


% ---------------------------
%Paramètres
% ---------------------------
epsilon=0.01*sqrt(2);
Nbpt_H=21;   %Maillage grossier
Tps_initial = 0;
Tps_final = 0.01*epsilon^2;
Nbpt_ref = 10001;  %Maillage fin
[lambda,Psi] = Solution_pb_spectral(Nbpt_ref);
delta_t = (epsilon^2)/(max(lambda,1)*2000);

Nombre_vect_propre = 15; %Nombre de vecteur propre que l'on calcul pour la reconstruction spectrale

Methode = 1; % Methode=0 --> Elements finis P1  /  Methode=1 --> Elements finis MsFEM
num_idee = 8 ; %Voir la description dans la fonction Khi_i() pour les détails des numéros (seulement dans le cas où Methode=1)

condition_initiale_artificielle = 1; %Si on construit la condition initiale artificiellement à l'aide des modes propres
Nombre_vect_propre_cond_initiale = 10;

affichage_erreurs = 1;
%Affichages des erreurs entre solution approchee et solution de référence :
% 0 Si pas d'affichage
% 1 pour affichage

trace_solution_Tf = 1;
%Tracé des solutions exactes et approchées rescalées à l'instant final :
% 0 Si pas de tracé
% 1 pour le tracé solution exacte et solution MsFEM

% ---------------------------


h_ref = 1/(Nbpt_ref -1);
alpha = 1/delta_t;
N_t = round((Tps_final-Tps_initial)/delta_t); % le nombre d'iterations necessaires
H = 1/(Nbpt_H -1);

% Maillage
% ----------------------
XX_ref = zeros(Nbpt_ref,1);     %vecteur des coordonnées des noeuds du maillage de pas H
for i=1:Nbpt_ref
    XX_ref(i)=(i-1)*h_ref;
end
XX_H = zeros(Nbpt_H,1);     %vecteur des coordonnées des noeuds du maillage de pas H
for i=1:Nbpt_H
    XX_H(i)=(i-1)*H;
end

% ----------------------------------------------
% Condition initiale
% ------------------
UU_0 = condition_initiale(XX_ref);
%UU_H = UU_0;

% calcul des matrices EF pour fonction de reference
% -------------------------------------------------
MM_ref=matM_ref(Nbpt_ref);
MM_sigma_ref=matM_sigma_ref(Nbpt_ref,epsilon);
KK_ref=matK_ref(Nbpt_ref,epsilon);
KK_ref_A_constant = matK_unitaire(Nbpt_ref,h_ref);

%AA_ref = alpha*MM_ref + (1/epsilon^2)*MM_sigma_ref + KK_ref;
AA_ref = alpha*MM_ref + KK_ref ;
AA_ref(1,1)=1; AA_ref(Nbpt_ref,Nbpt_ref)=1;

AA_stationnaire_ref = MM_sigma_ref + KK_ref*(epsilon^2);
AA_stationnaire_ref(1,1)=1; AA_stationnaire_ref(Nbpt_ref,Nbpt_ref)=1;

[EV_ref,DV_ref] = eigs(AA_stationnaire_ref,MM_ref,Nombre_vect_propre,'smallestabs');

if condition_initiale_artificielle
    UU_0 = zeros(Nbpt_ref,1);
    for k = Nombre_vect_propre_cond_initiale
        %UU_0 = UU_0 + 1/(sqrt(Nombre_vect_propre_cond_initiale)) .* EV_ref(:,k);
        UU_0 = UU_0 + EV_ref(:,k);
    end
end
UU_ref = UU_0;
% calcul des matrices EF pour solution approchee
% ----------------------------------------------
if Methode==0

    MM=matM_ref(Nbpt_H);
    MM_sigma=matM_sigma_H_P1(Nbpt_H,Nbpt_ref,epsilon);
    KK=matK_H_P1(Nbpt_H,Nbpt_ref,epsilon);

elseif Methode==1

    Nbpt_h = (Nbpt_ref - 1)/(Nbpt_H-1) +1;
    KK_ref_A_constant_loc = matK_unitaire(Nbpt_h-1,h_ref);

    Khi_global = zeros(Nbpt_h,Nbpt_H-2,2);
    for l=1:(Nbpt_H-2)
        [fct_forme_gauche , fct_forme_droite] = Khi_i(Nbpt_H,Nbpt_h,epsilon,l,num_idee);
    
        Khi_global(:,l,1)= fct_forme_gauche ;
        Khi_global(:,l,2)= fct_forme_droite ;
    end
    MM=matM(Khi_global,Nbpt_H,Nbpt_h);
    MM_sigma=matM_sigma(Khi_global,Nbpt_H,Nbpt_h,epsilon);
    KK=matK(Khi_global,Nbpt_H,Nbpt_h,epsilon);
    
    
    AA = MM_sigma + KK*(epsilon^2) ;
    AA(1,1)=1; AA(Nbpt_H,Nbpt_H)=1;
    
    [EV,DV] = eigs(AA,MM,Nombre_vect_propre,'smallestabs');
    %[EV,DV] = eigs(AA,MM,1,'smallestabs');
    
    Liste_UU_H_reconstruit_k = zeros(Nbpt_ref,Nombre_vect_propre);
    for k = 1:Nombre_vect_propre
        UU_H_reconstruit_k = zeros(Nbpt_ref,1);
        for i=2:(Nbpt_H-2)        %Ajout de la composante i du vecteur UU
            for j=1:(Nbpt_h-1)
                UU_H_reconstruit_k((i-1)*(Nbpt_h-1)+j) = UU_H_reconstruit_k((i-1)*(Nbpt_h-1)+j) + EV(i,k)*Khi_global(j,i-1,2) + EV(i+1,k)*Khi_global(j,i,1);
            end
        end
        %Premier intervalle
        i = 1; % On trace la solution entre 0 et H
        for j=1:(Nbpt_h-1)
            UU_H_reconstruit_k((i-1)*(Nbpt_h-1)+j) = UU_H_reconstruit_k((i-1)*(Nbpt_h-1)+j) + EV(i+1,k)*Khi_global(j,i,1);
        end
        %Dernier intervalle
        i = Nbpt_H-1; % On trace la solution entre 1-H et 1
        for j=1:Nbpt_h
            UU_H_reconstruit_k((i-1)*(Nbpt_h-1)+j) = UU_H_reconstruit_k((i-1)*(Nbpt_h-1)+j) + EV(i,k)*Khi_global(j,i-1,2);
        end
    Norme_L2 = sqrt(sum(UU_H_reconstruit_k.*UU_H_reconstruit_k)*h_ref);
    Liste_UU_H_reconstruit_k(:,k) = UU_H_reconstruit_k/Norme_L2;
    end

end

Produit_scalaireL2_ref = zeros(Nombre_vect_propre,1);
Produit_scalaireL2 = zeros(Nombre_vect_propre,1);
UU_H_reconstruit = zeros(Nbpt_ref,1);
Norme_L2_U0 = sqrt(sum(UU_0.*UU_0)*h_ref);
for k = 1:Nombre_vect_propre
    Produit_scalaireL2(k) = sum(UU_0.*Liste_UU_H_reconstruit_k(:,k))*h_ref;
    Produit_scalaireL2_ref(k) = abs(sum(UU_0.*EV_ref(:,k))*h_ref)/Norme_L2_U0;
    UU_H_reconstruit = UU_H_reconstruit + Produit_scalaireL2(k)*exp(-DV(k,k)*N_t*delta_t/epsilon^2).*Liste_UU_H_reconstruit_k(:,k);
end


Erreur_L2 = 0;
List_Erreur_L2_t=[];
Erreur_H1 = 0;
Erreur_H1_t = 0;
List_Erreur_H1_t=[];
for k = 1:N_t

    % Code pour chaque itération
    if mod(k, floor(N_t/10)) == 0 % Si le modulo est nul, c'est qu'on a fait un dixième de la boucle
        fprintf('%d%% terminé\n', round(100 * k / N_t)); % On affiche la progression en pourcentage
    end

    %Etape 0.5 : Partie réaction
    for i=1:Nbpt_ref
        xi = (i-1)*h_ref;
        UU_ref(i) = exp(-delta_t*Sigma(xi,epsilon)/(2*epsilon^2))*UU_ref(i); %On approxime l'exp par la valeur au noeud
    end

    % inversion
    % ----------
    LL_k_ref = alpha*MM_ref*UU_ref;

    UU_ref = AA_ref\LL_k_ref;

    %Etape 1.5 : Partie réaction
    for i=1:Nbpt_ref
        xi = (i-1)*h_ref;
        UU_ref(i) = exp(-delta_t*Sigma(xi,epsilon)/(2*epsilon^2))*UU_ref(i); %On approxime l'exp par la valeur au noeud
    end
end

% -------------------------------------------------------------------
% Calcul de toutes les erreurs
for i=1:(Nbpt_H-1)        %Ajout de la composante i du vecteur UU
    Vecteur_propre_ref_local = UU_ref((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h-1);
    solution_intervalle_i = UU_H_reconstruit((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h-1);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant_loc*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1_t = Erreur_H1_t + Erreur_H1_locale;
end

Erreur_H1_t = Erreur_H1_t/((UU_ref)'*KK_ref_A_constant*(UU_ref));
Erreur_L2_t = ((UU_H_reconstruit-UU_ref)'*MM_ref*(UU_H_reconstruit-UU_ref))/(UU_ref'*MM_ref*UU_ref);

if affichage_erreurs==1

    %Erreurs pour u_eps
    Erreur_L2_relative_instant_Tf = sqrt(Erreur_L2_t) %#ok<NOPTS>
    Erreur_H1_relative_instant_Tf = sqrt(Erreur_H1_t) %#ok<NOPTS>

end

% -------------------------------------------------------------------
% visualisation
if trace_solution_Tf ==1
    R = exp(Tps_final*lambda/epsilon^2);
    plot(XX_ref,R*UU_ref,'o')
    hold on
    plot(XX_ref,R*UU_H_reconstruit,LineWidth=2);
    axis([0,1, 0 , 1.2*R*max(UU_ref)]);
    legend('solution exacte','solution approchée')
    title(sprintf("T_f=%g / Dt=%g / eps=%g / H=%g",Tps_final,delta_t,epsilon,H))
end
%%%