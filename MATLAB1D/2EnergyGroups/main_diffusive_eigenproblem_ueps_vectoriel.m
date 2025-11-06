% =====================================================
%
% principal_reaction_diffusion pour le probleme en u_eps multigroup;
%
% une routine qui calcul pour un pas H donné l'erreur entre la solution approchée par la méthode de résolution de notre choix
% et la solution exacte, calculée avec les éléments P1 classiques de pas de
% maillage très petit

% ---------------------------
%Paramètres
% ---------------------------
epsilon=0.0572519083;
Nbpt_H=11;   %Maillage grossier
Nbpt_ref = 10001;  %Maillage fin

num_idee = 12 ; %Voir la description dans la fonction Khi_i() pour les détails des numéros
%L'idee MsFEM Qui fonctionne in fine est celle notée Numero 11

affichage_erreurs = 1;
%Affichages des erreurs entre solution approchee et solution de référence :
% 0 Si pas d'affichage
% 1 pour affichage

trace_solution_Tf = 1;
%Tracé des solutions exactes et approchées :
% 0 Si pas de tracé
% 1 pour le tracé solution exacte et solution MsFEM

% ---------------------------


h_ref = 1/(Nbpt_ref -1);
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
tic;

[val_propre_ref,Vecteur_propre_ref_1,Vecteur_propre_ref_2] = fonction_propre_ref_ueps(Nbpt_ref,epsilon);
normL2carre = sqrt(norm(Vecteur_propre_ref_1,2)^2/(Nbpt_ref-1) + norm(Vecteur_propre_ref_2,2)^2/(Nbpt_ref-1)); % Rescale pour une norme L^2 = 1
Vecteur_propre_ref_1 = Vecteur_propre_ref_1./normL2carre;
Vecteur_propre_ref_2 = Vecteur_propre_ref_2./normL2carre;

tempsEcoule = toc;
disp(['Le temps de calcul de la solution de reference est de ', num2str(tempsEcoule), ' secondes.']);
% calcul des matrices EF pour solution approchee
% ----------------------------------------------
tic;

Nbpt_h = (Nbpt_ref - 1)/(Nbpt_H-1) +1;
h = 1/(Nbpt_h -1);
Khi_global_1_cas1 = zeros(Nbpt_h,Nbpt_H-2,2);
Khi_global_1_cas2 = zeros(Nbpt_h,Nbpt_H-2,2);
Khi_global_2_cas1 = zeros(Nbpt_h,Nbpt_H-2,2);
Khi_global_2_cas2 = zeros(Nbpt_h,Nbpt_H-2,2);
for l=1:(Nbpt_H-2)
    [ UU_g_cas1_1, UU_g_cas1_2 , UU_d_cas1_1, UU_d_cas1_2, UU_g_cas2_1, UU_g_cas2_2 , UU_d_cas2_1, UU_d_cas2_2] = Khi_i_vectoriel(Nbpt_H,Nbpt_h,epsilon,l,num_idee);

    Khi_global_1_cas1(:,l,1)= UU_g_cas1_1 ;
    Khi_global_1_cas1(:,l,2)= UU_d_cas1_1 ;

    Khi_global_2_cas1(:,l,1)= UU_g_cas1_2 ;
    Khi_global_2_cas1(:,l,2)= UU_d_cas1_2 ;

    Khi_global_1_cas2(:,l,1)= UU_g_cas2_1 ;
    Khi_global_1_cas2(:,l,2)= UU_d_cas2_1 ;

    Khi_global_2_cas2(:,l,1)= UU_g_cas2_2 ;
    Khi_global_2_cas2(:,l,2)= UU_d_cas2_2 ;
end
%
% [~,Psi_1_adjoint,Psi_2_adjoint] = Solution_pb_spectral_adjoint(Nbpt_ref);
% [~,Psi_1,Psi_2] = Solution_pb_spectral(Nbpt_ref);
% Symmetry = check_symmetry(Psi_1,Psi_2,Psi_1_adjoint,Psi_2_adjoint);

Vect_ones = ones(Nbpt_ref,1);
MM_11=matM(Khi_global_1_cas1,Vect_ones , Vect_ones,Nbpt_H,Nbpt_h, epsilon) + matM(Khi_global_2_cas1,Vect_ones , Vect_ones,Nbpt_H,Nbpt_h, epsilon);
MM_22=matM(Khi_global_1_cas2,Vect_ones , Vect_ones,Nbpt_H,Nbpt_h, epsilon) + matM(Khi_global_2_cas2,Vect_ones , Vect_ones,Nbpt_H,Nbpt_h, epsilon);
MM_21=matM_vectoriel(Khi_global_1_cas2, Khi_global_1_cas1,Nbpt_H,Nbpt_h) + matM_vectoriel(Khi_global_2_cas2, Khi_global_2_cas1,Nbpt_H,Nbpt_h);

KK_11=matK(Khi_global_1_cas1,Vect_ones , Vect_ones,Nbpt_H,Nbpt_h,epsilon, 1) + matK(Khi_global_2_cas1,Vect_ones , Vect_ones,Nbpt_H,Nbpt_h,epsilon, 2);
KK_22=matK(Khi_global_1_cas2,Vect_ones , Vect_ones,Nbpt_H,Nbpt_h,epsilon, 1) + matK(Khi_global_2_cas2,Vect_ones , Vect_ones,Nbpt_H,Nbpt_h,epsilon, 2);

KK_21=matK_vectoriel(Khi_global_1_cas2, Khi_global_1_cas1,Vect_ones , Vect_ones, Nbpt_H, Nbpt_h, epsilon, 1) ...
    + matK_vectoriel(Khi_global_2_cas2, Khi_global_2_cas1,Vect_ones , Vect_ones, Nbpt_H, Nbpt_h, epsilon, 2);

MM_sigma_11=matM_sigma(Khi_global_1_cas1, Khi_global_1_cas1,Nbpt_H,Nbpt_h,epsilon,1) ...
    + matM_sigma(Khi_global_1_cas1, Khi_global_2_cas1,Nbpt_H,Nbpt_h,epsilon,2) ...
    + matM_sigma(Khi_global_2_cas1, Khi_global_1_cas1,Nbpt_H,Nbpt_h,epsilon,3) ...
    + matM_sigma(Khi_global_2_cas1, Khi_global_2_cas1,Nbpt_H,Nbpt_h,epsilon,4);

MM_sigma_22=matM_sigma(Khi_global_1_cas2, Khi_global_1_cas2,Nbpt_H,Nbpt_h,epsilon,1) ...
    + matM_sigma(Khi_global_1_cas2, Khi_global_2_cas2,Nbpt_H,Nbpt_h,epsilon,2) ...
    + matM_sigma(Khi_global_2_cas2, Khi_global_1_cas2,Nbpt_H,Nbpt_h,epsilon,3) ...
    + matM_sigma(Khi_global_2_cas2, Khi_global_2_cas2,Nbpt_H,Nbpt_h,epsilon,4);

MM_sigma_21=matM_sigma(Khi_global_1_cas2, Khi_global_1_cas1,Nbpt_H,Nbpt_h,epsilon,1) ...
    + matM_sigma(Khi_global_1_cas2, Khi_global_2_cas1,Nbpt_H,Nbpt_h,epsilon,2) ...
    + matM_sigma(Khi_global_2_cas2, Khi_global_1_cas1,Nbpt_H,Nbpt_h,epsilon,3) ...
    + matM_sigma(Khi_global_2_cas2, Khi_global_2_cas1,Nbpt_H,Nbpt_h,epsilon,4);

MM_sigma_12=matM_sigma(Khi_global_1_cas1, Khi_global_1_cas2,Nbpt_H,Nbpt_h,epsilon,1) ...
    + matM_sigma(Khi_global_1_cas1, Khi_global_2_cas2,Nbpt_H,Nbpt_h,epsilon,2) ...
    + matM_sigma(Khi_global_2_cas1, Khi_global_1_cas2,Nbpt_H,Nbpt_h,epsilon,3) ...
    + matM_sigma(Khi_global_2_cas1, Khi_global_2_cas2,Nbpt_H,Nbpt_h,epsilon,4);

%Assemblage de toutes les matrices

Zeros = sparse(Nbpt_H,Nbpt_H); % matrice de masse

KK = [KK_11 , KK_21' ; KK_21 ,KK_22 ];
MM_sigma = [MM_sigma_11 , MM_sigma_12 ; MM_sigma_21 , MM_sigma_22];
MM = [MM_11 , MM_21' ; MM_21 ,MM_22 ];

AA = KK*(epsilon^2) + MM_sigma;
AA(1,1)=1; AA(2*Nbpt_H,2*Nbpt_H)=1;
AA(Nbpt_H,Nbpt_H)=1 ; AA(Nbpt_H+1,Nbpt_H+1)=1;

%[EV,DV] = eigs(AA,MM,5,'smallestabs');
[EV,DV] = eigs(AA,MM,1,'smallestabs');


[Valeur_propre_min_H,arg_min] = min(diag(DV));
Vecteur_propre_H = EV(:,arg_min);

normL2 = sqrt(Vecteur_propre_H'*MM*Vecteur_propre_H);
Vecteur_propre_H = Vecteur_propre_H/normL2;

Vecteur_propre_H_1 = Vecteur_propre_H(1:Nbpt_H);
Vecteur_propre_H_2 = Vecteur_propre_H((Nbpt_H+1):2*Nbpt_H);

%On a un VP de norme L^2 unitaire. On choisis le VP positif.
if (Vecteur_propre_H_1(fix(Nbpt_H/2))<0)
    Vecteur_propre_H_1 = -Vecteur_propre_H_1;
end
if (Vecteur_propre_H_2(fix(Nbpt_H/2))<0)
    Vecteur_propre_H_2 = -Vecteur_propre_H_2;
end

% % Rescale pour une norme L^2 = 1
% MM_1_constant=matM(Khi_global_1, ones(Nbpt_ref,1) , ones(Nbpt_ref,1),Nbpt_H,Nbpt_h, epsilon);
% MM_2_constant=matM(Khi_global_2, ones(Nbpt_ref,1) , ones(Nbpt_ref,1),Nbpt_H,Nbpt_h, epsilon);
% 
% Norme_L2_1 = sqrt(Vecteur_propre_H_1'*MM_1_constant*Vecteur_propre_H_1);
% Vecteur_propre_H_1 = Vecteur_propre_H_1/Norme_L2_1;
% 
% Norme_L2_2 = sqrt(Vecteur_propre_H_2'*MM_2_constant*Vecteur_propre_H_2);
% Vecteur_propre_H_2 = Vecteur_propre_H_2/Norme_L2_2;

tempsEcoule = toc;
disp(['Le temps de calcul de la solution approchee est de ', num2str(tempsEcoule), ' secondes.']);

if trace_solution_Tf ==1
    figure;
    plot(XX_ref,Vecteur_propre_ref_1,'o')
    hold on ;

    for i=2:(Nbpt_H-2)
        % On trace la solution entre (i-1)H et iH
        solution_intervalle_i = zeros(Nbpt_ref,1);
        solution_intervalle_i = solution_intervalle_i*NaN;
        for j=1:Nbpt_h
            solution_intervalle_i((i-1)*(Nbpt_h-1)+j) = Vecteur_propre_H_1(i)*Khi_global_1_cas1(j,i-1,2) + Vecteur_propre_H_1(i+1)*Khi_global_1_cas1(j,i,1)...
                                     + Vecteur_propre_H_2(i)*Khi_global_1_cas2(j,i-1,2) + Vecteur_propre_H_2(i+1)*Khi_global_1_cas2(j,i,1);
        end
        plot(XX_ref,solution_intervalle_i,LineWidth=2,color='r')
        hold on;
    end

    %Premier intervalle
    solution_intervalle_i = zeros(Nbpt_ref,1);
    solution_intervalle_i = solution_intervalle_i*NaN;
    i = 1; % On trace la solution entre 0 et H
    for j=1:Nbpt_h
        solution_intervalle_i((i-1)*(Nbpt_h-1)+j) = Vecteur_propre_H_1(i+1)*Khi_global_1_cas1(j,i,1) + Vecteur_propre_H_2(i+1)*Khi_global_1_cas2(j,i,1);
    end
    plot(XX_ref,solution_intervalle_i,LineWidth=2,color='r')

    %Dernier intervalle
    solution_intervalle_i = zeros(Nbpt_ref,1);
    solution_intervalle_i = solution_intervalle_i*NaN;
    i = Nbpt_H-1; % On trace la solution entre 1-H et 1
    for j=1:Nbpt_h
        solution_intervalle_i((i-1)*(Nbpt_h-1)+j) = Vecteur_propre_H_1(i)*Khi_global_1_cas1(j,i-1,2) + Vecteur_propre_H_2(i)*Khi_global_1_cas2(j,i-1,2);
    end
    plot(XX_ref,solution_intervalle_i,LineWidth=2,color='r')

    axis([0,1, 0 , 1.2*max(Vecteur_propre_ref_1)]);
    legend('solution exacte Energie 1','solution approchée Energie 1')
    title(sprintf("eps=%g / H=%g",epsilon,H))
    set(gca,'FontSize',30)
    hold off;
end

if trace_solution_Tf ==1
    figure;
    plot(XX_ref,Vecteur_propre_ref_2,'o')
    hold on ;

    for i=2:(Nbpt_H-2)
        % On trace la solution entre (i-1)H et iH
        solution_intervalle_i = zeros(Nbpt_ref,1);
        solution_intervalle_i = solution_intervalle_i*NaN;
        for j=1:Nbpt_h
            solution_intervalle_i((i-1)*(Nbpt_h-1)+j) = Vecteur_propre_H_1(i)*Khi_global_2_cas1(j,i-1,2) + Vecteur_propre_H_1(i+1)*Khi_global_2_cas1(j,i,1)...
                                     + Vecteur_propre_H_2(i)*Khi_global_2_cas2(j,i-1,2) + Vecteur_propre_H_2(i+1)*Khi_global_2_cas2(j,i,1);
        end
        plot(XX_ref,solution_intervalle_i,LineWidth=2,color='r')
        hold on;
    end

    %Premier intervalle
    solution_intervalle_i = zeros(Nbpt_ref,1);
    solution_intervalle_i = solution_intervalle_i*NaN;
    i = 1; % On trace la solution entre 0 et H
    for j=1:Nbpt_h
        solution_intervalle_i((i-1)*(Nbpt_h-1)+j) = Vecteur_propre_H_1(i+1)*Khi_global_2_cas1(j,i,1) + Vecteur_propre_H_2(i+1)*Khi_global_2_cas2(j,i,1);
    end
    plot(XX_ref,solution_intervalle_i,LineWidth=2,color='r')

    %Dernier intervalle
    solution_intervalle_i = zeros(Nbpt_ref,1);
    solution_intervalle_i = solution_intervalle_i*NaN;
    i = Nbpt_H-1; % On trace la solution entre 1-H et 1
    for j=1:Nbpt_h
        solution_intervalle_i((i-1)*(Nbpt_h-1)+j) = Vecteur_propre_H_1(i)*Khi_global_2_cas1(j,i-1,2) + Vecteur_propre_H_2(i)*Khi_global_2_cas2(j,i-1,2);
    end
    plot(XX_ref,solution_intervalle_i,LineWidth=2,color='r')

    axis([0,1, 0 , 1.2*max(Vecteur_propre_ref_2)]);
    legend('solution exacte Energie 2','solution approchée Energie 2')
    title(sprintf("eps=%g / H=%g",epsilon,H))
    set(gca,'FontSize',30)
    hold off;
end
% -------------------------------------------------------------------
% Calcul de toutes les erreurs
if affichage_erreurs==1
    Erreur_H1_1 = 0;
    KK_ref_A_constant = matK_unitaire(Nbpt_h,h_ref);
    Norme_H1_vect_propre_ref_1 = 0;

    for i=2:(Nbpt_H-2)
        Vecteur_propre_ref_local = Vecteur_propre_ref_1((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
        solution_intervalle_i = zeros(Nbpt_h,1);
        for j=1:Nbpt_h
            solution_intervalle_i(j) = Vecteur_propre_H_1(i)*Khi_global_1_cas1(j,i-1,2) + Vecteur_propre_H_1(i+1)*Khi_global_1_cas1(j,i,1)...
                                     + Vecteur_propre_H_2(i)*Khi_global_1_cas2(j,i-1,2) + Vecteur_propre_H_2(i+1)*Khi_global_1_cas2(j,i,1);
        end
        Norme_H1_vect_propre_ref_1 = Norme_H1_vect_propre_ref_1 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
        Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
        Erreur_H1_1 = Erreur_H1_1 + Erreur_H1_locale;
    end

    %Premier intervalle
    i = 1; % On trace la solution entre 0 et H
    Vecteur_propre_ref_local = Vecteur_propre_ref_1((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h
        solution_intervalle_i(j) = Vecteur_propre_H_1(i+1)*Khi_global_1_cas1(j,i,1) + Vecteur_propre_H_2(i+1)*Khi_global_1_cas2(j,i,1);
    end
    Norme_H1_vect_propre_ref_1 = Norme_H1_vect_propre_ref_1 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1_1 = Erreur_H1_1 + Erreur_H1_locale;

    %Dernier intervalle
    i = Nbpt_H-1; % On trace la solution entre 1-H et 1
    Vecteur_propre_ref_local = Vecteur_propre_ref_1((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h
        solution_intervalle_i(j) = Vecteur_propre_H_1(i)*Khi_global_1_cas1(j,i-1,2) + Vecteur_propre_H_2(i)*Khi_global_1_cas2(j,i-1,2);
    end
    Norme_H1_vect_propre_ref_1 = Norme_H1_vect_propre_ref_1 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1_1 = Erreur_H1_1 + Erreur_H1_locale;

    Erreur_H1_relative_1 = Erreur_H1_1/Norme_H1_vect_propre_ref_1;
    Erreur_H1_relative_1_tmp = sqrt(Erreur_H1_relative_1) %#ok<NOPTS>

    Erreur_H1_2 = 0;
    KK_ref_A_constant = matK_unitaire(Nbpt_h,h_ref);
    Norme_H1_vect_propre_ref_2 = 0;

    for i=2:(Nbpt_H-2)
        Vecteur_propre_ref_local = Vecteur_propre_ref_2((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
        solution_intervalle_i = zeros(Nbpt_h,1);
        for j=1:Nbpt_h
            solution_intervalle_i(j) = Vecteur_propre_H_1(i)*Khi_global_2_cas1(j,i-1,2) + Vecteur_propre_H_1(i+1)*Khi_global_2_cas1(j,i,1)...
                                     + Vecteur_propre_H_2(i)*Khi_global_2_cas2(j,i-1,2) + Vecteur_propre_H_2(i+1)*Khi_global_2_cas2(j,i,1);
        end
        Norme_H1_vect_propre_ref_2 = Norme_H1_vect_propre_ref_2 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
        Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
        Erreur_H1_2 = Erreur_H1_2 + Erreur_H1_locale;
    end

    %Premier intervalle
    i = 1; % On trace la solution entre 0 et H
    Vecteur_propre_ref_local = Vecteur_propre_ref_2((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h
        solution_intervalle_i(j) = Vecteur_propre_H_1(i+1)*Khi_global_2_cas1(j,i,1) + Vecteur_propre_H_2(i+1)*Khi_global_2_cas2(j,i,1);
    end
    Norme_H1_vect_propre_ref_2 = Norme_H1_vect_propre_ref_2 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1_2 = Erreur_H1_2 + Erreur_H1_locale;

    %Dernier intervalle
    i = Nbpt_H-1; % On trace la solution entre 1-H et 1
    Vecteur_propre_ref_local = Vecteur_propre_ref_2((i-1)*(Nbpt_h-1)+1:(i-1)*(Nbpt_h-1)+Nbpt_h);
    solution_intervalle_i = zeros(Nbpt_h,1);
    for j=1:Nbpt_h
        solution_intervalle_i(j) = Vecteur_propre_H_1(i)*Khi_global_2_cas1(j,i-1,2) + Vecteur_propre_H_2(i)*Khi_global_2_cas2(j,i-1,2);
    end
    Norme_H1_vect_propre_ref_2 = Norme_H1_vect_propre_ref_2 + (Vecteur_propre_ref_local'*KK_ref_A_constant*Vecteur_propre_ref_local);
    Erreur_H1_locale = ((solution_intervalle_i-Vecteur_propre_ref_local)'*KK_ref_A_constant*(solution_intervalle_i-Vecteur_propre_ref_local));
    Erreur_H1_2 = Erreur_H1_2 + Erreur_H1_locale;

    Erreur_H1_relative_2 = Erreur_H1_2/Norme_H1_vect_propre_ref_2;
    Erreur_H1_relative_2_tmp = sqrt(Erreur_H1_relative_2) %#ok<NOPTS>



    Erreur_H1_relative = (1/sqrt(2))*sqrt(Erreur_H1_relative_1 + Erreur_H1_relative_2) %#ok<NOPTS>
    
    Erreur_VP = abs(val_propre_ref - Valeur_propre_min_H)/val_propre_ref
end
%%%
