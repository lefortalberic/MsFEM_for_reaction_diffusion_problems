%Code qui compare et trace les erreurs pour differentes methodes, en fonction de epsilon
% ---------------------------
%Paramètres
% ---------------------------
Nbpt_H=9;
H=1/(Nbpt_H -1);
Nbpt_H_list = 9;
epsilon = [0.1250,0.10714, 0.0833333333,  0.075, 0.068181, 0.0576923, 0.05,0.04411764, 0.04166, 0.039473,0.03 ,0.0241935,0.01415094,0.0105633];
Nbpt_ref = 2001;  %Maillage fin

h_ref = 1/(Nbpt_ref -1);

List_erreur_H1_P1 = [];
List_erreur_H1_triche = [];
List_erreur_H1_0_3H_periodique = [];
List_erreur_H1_0_3H_periodique_filtre = [];

%Solutions MsFEM
for eps=epsilon 

    %Solution de référence
    fprintf('\n Calcul de la solution de reference \n');
    tic;
    [~,Vecteur_propre_ref] = fonction_propre_eps_Dirichlet(2,Nbpt_ref,eps,1);
    tempsEcoule = toc;
    disp(['Le temps d''exécution de la solution de reference est de ', num2str(tempsEcoule), ' secondes.']);
    
    erreur_H1_P1 = Liste_Erreur_H1_relative(6,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref);
    erreur_H1_triche = Liste_Erreur_H1_relative(8,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref);
    erreur_H1_0_3H_periodique = Liste_Erreur_H1_relative(19,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref);
    erreur_H1_0_3H_periodique_filtre = Liste_Erreur_H1_relative(23,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref);

    List_erreur_H1_P1 = [List_erreur_H1_P1,erreur_H1_P1]; %#ok<*AGROW> 
    List_erreur_H1_triche = [List_erreur_H1_triche,erreur_H1_triche];
    List_erreur_H1_0_3H_periodique = [List_erreur_H1_0_3H_periodique,min(erreur_H1_0_3H_periodique,1)];
    List_erreur_H1_0_3H_periodique_filtre = [List_erreur_H1_0_3H_periodique_filtre,erreur_H1_0_3H_periodique_filtre];

end

% visualisation
% ------------------------------------------------------------------------------

h_list = 3*H./epsilon;
plot((h_list),(List_erreur_H1_P1),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot((h_list),(List_erreur_H1_triche),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot((h_list),(List_erreur_H1_0_3H_periodique),'LineWidth',4,'marker','x','markersize',28, 'Color', [1, 0.8, 0])
hold on;
plot((h_list),(List_erreur_H1_0_3H_periodique_filtre),'LineWidth',4,'marker','x','markersize',28, 'Color', [0.5, 0, 0.5])
hold on;

% Tracé des lignes verticales
ylims = get(gca,'YLim');
hold on;
plot([3*H/(0.03) 3*H/(0.03)],ylims,'--','Color','black','LineWidth',3);
text(3*H/(0.03),-0.03, ' \epsilon_{Min,2D}', 'HorizontalAlignment', 'center','FontSize', 28);
hold on;

legend('P1 method','Methode triche','MsFEM periodique 3H method','MsFEM 3H filtre','FontSize',25)

xlabel("3H/\epsilon",'FontSize',25)
ylabel("Err_{relative}",'FontSize',25)
set(gca,'FontSize',30)
hold off;
