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
Nbpt_H=11;   %Maillage grossier
Tps_initial = 0;
Tps_final = 6*epsilon^2;
Nbpt_ref = 10001;  %Maillage fin
[lambda,Psi] = Solution_pb_spectral(Nbpt_ref);
delta_t = (epsilon^2)/(max(lambda,1)*2000);

Methode = 1; % Methode=0 --> Elements finis P1  /  Methode=1 --> Elements finis MsFEM
num_idee = 8 ; %Voir la description dans la fonction Khi_i() pour les détails des numéros (seulement dans le cas où Methode=1)

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
UU_ref = condition_initiale(XX_ref);
UU_H = condition_initiale(XX_H);

% calcul des matrices EF pour fonction de reference
% -------------------------------------------------
MM_ref=matM_ref(Nbpt_ref);
%MM_sigma_ref=matM_sigma_ref(Nbpt_ref,epsilon);
KK_ref=matK_ref(Nbpt_ref,epsilon);
KK_ref_A_constant = matK_unitaire(Nbpt_ref,h_ref);

%AA_ref = alpha*MM_ref + (1/epsilon^2)*MM_sigma_ref + KK_ref ;
AA_ref = alpha*MM_ref + KK_ref ;
AA_ref(1,1)=1; AA_ref(Nbpt_ref,Nbpt_ref)=1;

% calcul des matrices EF pour solution approchee
% ----------------------------------------------
if Methode==0

    MM=matM_ref(Nbpt_H);
    MM_sigma=matM_sigma_H_P1(Nbpt_H,Nbpt_ref,epsilon);
    KK=matK_H_P1(Nbpt_H,Nbpt_ref,epsilon);

elseif Methode==1

    Nbpt_h = (Nbpt_ref - 1)/(Nbpt_H-1) +1;
    KK_ref_A_constant_loc = matK_unitaire(Nbpt_h,h_ref);
    Khi_global = zeros(Nbpt_h,Nbpt_H-2,2);
    for l=1:(Nbpt_H-2)
        [fct_forme_gauche , fct_forme_droite] = Khi_i(Nbpt_H,Nbpt_h,epsilon,l,num_idee);
    
        Khi_global(:,l,1)= fct_forme_gauche ;
        Khi_global(:,l,2)= fct_forme_droite ;
    end
    MM=matM(Khi_global,Nbpt_H,Nbpt_h);
    MM_sigma=matM_sigma(Khi_global,Nbpt_H,Nbpt_h,epsilon);
    KK=matK(Khi_global,Nbpt_H,Nbpt_h,epsilon);
    
%     for i=1:(Nbpt_H-2)
%         UU_H(i+1) = UU_H(i+1)/Khi_global(Nbpt_h,i); %On rescale la condition initiale par les fonctions de forme
%     end
end

AA = alpha*MM + (1/epsilon^2)*MM_sigma + KK ;
AA(1,1)=1; AA(Nbpt_H,Nbpt_H)=1;

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
    LL_k_H = alpha*MM*UU_H;

    UU_ref = AA_ref\LL_k_ref;
    UU_H = AA\LL_k_H;    %UU vaut ensuite UU_(t=k*delta_t)

    %Etape 1.5 : Partie réaction
    for i=1:Nbpt_ref
        xi = (i-1)*h_ref;
        UU_ref(i) = exp(-delta_t*Sigma(xi,epsilon)/(2*epsilon^2))*UU_ref(i); %On approxime l'exp par la valeur au noeud
    end

    % Projection sur la base des fonctions de formes (Pour P1 c'est juste
    % une interpolation)
    % -------------------------------------------------------------------
    if Methode==0
        UU_H_projetee = decompose(UU_H,Nbpt_ref);%On interpole UU_P1 (ie maillage grossier), sur un maillage plus fin
    elseif Methode==1
        solution_intervalle_i = zeros(Nbpt_h,1);
        UU_H_projetee = zeros(Nbpt_ref,1);
        for i=2:(Nbpt_H-2)        %Ajout de la composante i du vecteur UU
            Vecteur_propre_ref_local = UU_ref((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
            for j=1:(Nbpt_h-1)
                UU_H_projetee((i-1)*(Nbpt_h-1)+j) = UU_H_projetee((i-1)*(Nbpt_h-1)+j) + UU_H(i)*Khi_global(j,i-1,2) + UU_H(i+1)*Khi_global(j,i,1);
                solution_intervalle_i(j) = UU_H_projetee((i-1)*(Nbpt_h-1)+j);
            end
            j= Nbpt_h;
            solution_intervalle_i(j) = UU_H_projetee((i-1)*(Nbpt_h-1)+j) + UU_H(i)*Khi_global(j,i-1,2) + UU_H(i+1)*Khi_global(j,i,1);
            Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant_loc*(solution_intervalle_i-Vecteur_propre_ref_local));
            Erreur_H1_t = Erreur_H1_t + Erreur_H1_locale;
        end
        %Premier intervalle
        i = 1; % On trace la solution entre 0 et H
        Vecteur_propre_ref_local = UU_ref((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
        solution_intervalle_i = zeros(Nbpt_h,1);
        for j=1:(Nbpt_h-1)
            UU_H_projetee((i-1)*(Nbpt_h-1)+j) = UU_H_projetee((i-1)*(Nbpt_h-1)+j) + UU_H(i+1)*Khi_global(j,i,1);
            solution_intervalle_i(j) = UU_H_projetee((i-1)*(Nbpt_h-1)+j);
        end
        j= Nbpt_h;
        solution_intervalle_i(j) = UU_H(i+1)*Khi_global(j,i,1);
        Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant_loc*(solution_intervalle_i-Vecteur_propre_ref_local));
        Erreur_H1_t = Erreur_H1_t + Erreur_H1_locale;
        %Dernier intervalle
        i = Nbpt_H-1; % On trace la solution entre 1-H et 1
        Vecteur_propre_ref_local = UU_ref((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
        solution_intervalle_i = zeros(Nbpt_h,1);
        for j=1:Nbpt_h
            UU_H_projetee((i-1)*(Nbpt_h-1)+j) = UU_H_projetee((i-1)*(Nbpt_h-1)+j) + UU_H(i)*Khi_global(j,i-1,2);
            solution_intervalle_i(j) = UU_H_projetee((i-1)*(Nbpt_h-1)+j);
        end
        Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant_loc*(solution_intervalle_i-Vecteur_propre_ref_local));
        Erreur_H1_t = Erreur_H1_t + Erreur_H1_locale;

        Erreur_H1_t = Erreur_H1_t/((UU_ref)'*KK_ref_A_constant*(UU_ref));
    end
    % -------------------------------------------------------------------

    Erreur_L2_t = ((UU_H_projetee-UU_ref)'*MM_ref*(UU_H_projetee-UU_ref))/(UU_ref'*MM_ref*UU_ref);

    List_Erreur_L2_t = [List_Erreur_L2_t ; sqrt(Erreur_L2_t)]; %#ok<AGROW>
    List_Erreur_H1_t = [List_Erreur_H1_t ; sqrt(Erreur_H1_t)]; %#ok<AGROW>

    Erreur_L2 = Erreur_L2 + Erreur_L2_t;
    Erreur_H1 = Erreur_H1 + Erreur_H1_t;
    Erreur_H1_t = 0;
end

% -------------------------------------------------------------------
% Calcul de toutes les erreurs
if affichage_erreurs==1

    %Erreurs pour u_eps
    Erreur_L2_relative_moyenne_en_temps_sur_0_Tf = sqrt(Erreur_L2/N_t) %#ok<NOPTS>
    Erreur_L2_relative_instant_Tf = List_Erreur_L2_t(N_t) %#ok<NOPTS>
    Erreur_L2_relative_a_DeltaT = List_Erreur_L2_t(1) %#ok<NOPTS>

    Erreur_H1_relative_moyenne_en_temps_sur_0_Tf = sqrt(Erreur_H1/N_t) %#ok<NOPTS>
    Erreur_H1_relative_instant_Tf = List_Erreur_H1_t(N_t) %#ok<NOPTS>
    Erreur_H1_relative_a_DeltaT = List_Erreur_H1_t(1) %#ok<NOPTS>

    %Erreurs pour u_eps/Psi_eps
    VV_eps = zeros(Nbpt_ref,1);
    for i=1:Nbpt_ref
        xi = (i-1)*h_ref ;
        yi = rem(xi/epsilon,1) ; % y vaut x/epsilon modulo 1
        yk = fix(yi*(Nbpt_ref-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)
        VV_eps(i) = Psi(yk);
    end

    UU_H_projetee_sur_Psi_eps_exacte = UU_H_projetee./VV_eps ;
    UU_ref_sur_Psi_eps_exacte = UU_ref./VV_eps ;

    Erreur_L2_t_u_sur_Psi_exacte = ((UU_H_projetee_sur_Psi_eps_exacte-UU_ref_sur_Psi_eps_exacte)'*MM_ref*(UU_H_projetee_sur_Psi_eps_exacte-UU_ref_sur_Psi_eps_exacte))/(UU_ref_sur_Psi_eps_exacte'*MM_ref*UU_ref_sur_Psi_eps_exacte);
    Erreur_H1_t_u_sur_Psi_exacte = ((UU_H_projetee_sur_Psi_eps_exacte-UU_ref_sur_Psi_eps_exacte)'*KK_ref_A_constant*(UU_H_projetee_sur_Psi_eps_exacte-UU_ref_sur_Psi_eps_exacte))/(UU_ref_sur_Psi_eps_exacte'*KK_ref_A_constant*UU_ref_sur_Psi_eps_exacte);

    Erreur_L2_relative_instant_Tf_ueps_sur_Psieps_exacte = sqrt(Erreur_L2_t_u_sur_Psi_exacte) %#ok<NOPTS>

    Erreur_H1_relative_instant_Tf_ueps_sur_Psieps_exacte = sqrt(Erreur_H1_t_u_sur_Psi_exacte) %#ok<NOPTS>

    if (num_idee==9 && Methode==1)
        Nbpt_pour_mauvaise_approximation = 11 ; %Valeur à modifier en fonction de l'approximation plus ou moins mauvaise que l'on veut
        [~,Psi_mauvaise_approximation] = Solution_pb_spectral(Nbpt_pour_mauvaise_approximation);
        Psi_mauvaise_approximation = decompose(Psi_mauvaise_approximation,Nbpt_ref); %Interpollation de sur le maillage fin de l'approximation grossière de Psi
        Psi_eps_mauvaise_approximation = zeros(Nbpt_ref,1);

        for i=1:Nbpt_ref
            xi = (i-1)*h_ref ;
            yi = rem(xi/epsilon,1) ; % y vaut x/epsilon modulo 1
            yk = fix(yi*(Nbpt_ref-1)) +1; %partie entière de y*(Nbpt-1) (indice pour les matrices)

            Psi_eps_mauvaise_approximation(i) = Psi_mauvaise_approximation(yk);
            %t = yi*(Nbpt_ref-1) - fix(yi*(Nbpt_ref-1)) ;
            %Psi_eps_mauvaise_approximation(i) = Psi_mauvaise_approximation(yk)*(1-t) +Psi_mauvaise_approximation(yk+1)*t;
        end

        UU_H_projetee_sur_Psi_eps_approx = UU_H_projetee./Psi_eps_mauvaise_approximation ;
        UU_ref_sur_Psi_eps_approx = UU_ref./Psi_eps_mauvaise_approximation ;

        Erreur_L2_t_u_sur_Psi_approx = ((UU_H_projetee_sur_Psi_eps_approx-UU_ref_sur_Psi_eps_approx)'*MM_ref*(UU_H_projetee_sur_Psi_eps_approx-UU_ref_sur_Psi_eps_approx))/(UU_ref_sur_Psi_eps_approx'*MM_ref*UU_ref_sur_Psi_eps_approx);
        Erreur_H1_t_u_sur_Psi_approx = ((UU_H_projetee_sur_Psi_eps_approx-UU_ref_sur_Psi_eps_approx)'*KK_ref_A_constant*(UU_H_projetee_sur_Psi_eps_approx-UU_ref_sur_Psi_eps_approx))/(UU_ref_sur_Psi_eps_approx'*KK_ref_A_constant*UU_ref_sur_Psi_eps_approx);

        Erreur_L2_relative_instant_Tf_ueps_sur_Psieps_approx = sqrt(Erreur_L2_t_u_sur_Psi_approx) %#ok<NOPTS>

        Erreur_H1_relative_instant_Tf_ueps_sur_Psieps_approx = sqrt(Erreur_H1_t_u_sur_Psi_approx) %#ok<NOPTS>

        Erreur_L2_t_u_sur_Psi_semi_approx = ((UU_H_projetee_sur_Psi_eps_approx-UU_ref_sur_Psi_eps_exacte)'*MM_ref*(UU_H_projetee_sur_Psi_eps_approx-UU_ref_sur_Psi_eps_exacte))/(UU_ref_sur_Psi_eps_exacte'*MM_ref*UU_ref_sur_Psi_eps_exacte);
        Erreur_H1_t_u_sur_Psi_semi_approx = ((UU_H_projetee_sur_Psi_eps_approx-UU_ref_sur_Psi_eps_exacte)'*KK_ref_A_constant*(UU_H_projetee_sur_Psi_eps_approx-UU_ref_sur_Psi_eps_exacte))/(UU_ref_sur_Psi_eps_exacte'*KK_ref_A_constant*UU_ref_sur_Psi_eps_exacte);

        Erreur_L2_relative_instant_Tf_ueps_sur_Psieps_semi_approx = sqrt(Erreur_L2_t_u_sur_Psi_semi_approx) %#ok<NOPTS>

        Erreur_H1_relative_instant_Tf_ueps_sur_Psieps_semi_approx = sqrt(Erreur_H1_t_u_sur_Psi_semi_approx) %#ok<NOPTS>
    end
end

% -------------------------------------------------------------------
% visualisation
if trace_solution_Tf ==1
    R = exp(Tps_final*lambda/epsilon^2);
    plot(XX_ref,R*UU_ref,'o')
    hold on
    plot(XX_ref,R*UU_H_projetee,LineWidth=2);
    axis([0,1, 0 , 1.2*R*max(UU_ref)]);
    legend('solution exacte','solution approchée')
    title(sprintf("T_f=%g / Dt=%g / eps=%g / H=%g",Tps_final,delta_t,epsilon,H))
end
%%%