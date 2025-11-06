% ---------------------------
%Paramètres
% ---------------------------

%Nbpt_H_list=[3,5,6,11,21,26,51,101,126,201,251,501,1001];
Nbpt_H_list=[9];
H=1/(Nbpt_H_list -1);
%Nbpt_H_list=[5,11,26,51,101,126,201,251];
%Nbpt_H_list = 251;
%epsilon = [0.1250,0.10714, 0.0833333333,  0.075, 0.068181, 0.0576923, 0.05,0.04411764, 0.04166, 0.039473,0.03 ,0.0241935,0.01415094,0.0105633 ,0.00572519083,0.00299281723];
epsilon = [0.1250,0.10714, 0.0833333333,  0.075, 0.068181, 0.0576923, 0.05,0.04411764, 0.04166, 0.039473,0.03 ,0.0241935,0.01415094,0.0105633];
epsilon = [0.1250,0.10714, 0.0833333333,  0.075, 0.068181, 0.0576923, 0.05,0.04411764, 0.04166, 0.039473,0.03 ,0.0241935];

Nbpt_ref = 3001;  %Maillage fin

h_ref = 1/(Nbpt_ref -1);

List_erreur_H1_P1 = [];
List_erreur_H1_0_3H = [];
List_erreur_H1_0_3H_periodique = [];
List_erreur_H1_0_9H_periodique = [];
List_erreur_H1_0_9H = [];

%Solutions MsFEM
for eps=epsilon 

    %Solution de référence
    fprintf('\n Calcul de la solution de reference \n');
    tic;
    [~,Vecteur_propre_ref] = fonction_propre_eps_Dirichlet(2,Nbpt_ref,eps,1);
    tempsEcoule = toc;
    disp(['Le temps d''exécution de la solution de reference est de ', num2str(tempsEcoule), ' secondes.']);
    
    erreur_H1_P1 = Liste_Erreur_H1_relative(6,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref);
    erreur_H1_0_3H = Liste_Erreur_H1_relative(11,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref);
    erreur_H1_0_3H_periodique = Liste_Erreur_H1_relative(19,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref);
    erreur_H1_0_9H_periodique = Liste_Erreur_H1_relative(21,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref);
    erreur_H1_0_9H = Liste_Erreur_H1_relative(16,Nbpt_H_list,eps,Nbpt_ref,Vecteur_propre_ref);

    List_erreur_H1_P1 = [List_erreur_H1_P1,erreur_H1_P1]; %#ok<*AGROW> 
    List_erreur_H1_0_3H = [List_erreur_H1_0_3H,erreur_H1_0_3H];
    List_erreur_H1_0_3H_periodique = [List_erreur_H1_0_3H_periodique,min(erreur_H1_0_3H_periodique,1)];
    List_erreur_H1_0_9H_periodique = [List_erreur_H1_0_9H_periodique,min(erreur_H1_0_9H_periodique,1)];
    List_erreur_H1_0_9H = [List_erreur_H1_0_9H,erreur_H1_0_9H];

end

% visualisation
% ------------------------------------------------------------------------------

h_list = 3*H./epsilon;
plot((h_list),(List_erreur_H1_P1),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot((h_list),(List_erreur_H1_0_3H),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot((h_list),(List_erreur_H1_0_9H),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot((h_list),(List_erreur_H1_0_3H_periodique),'LineWidth',4,'marker','x','markersize',28)
hold on;
plot((h_list),(List_erreur_H1_0_9H_periodique),'LineWidth',4,'marker','x','markersize',28)
hold on;
%plot(-log10(h_list),log10(h_list.^1)-0.15,'LineWidth',4)
hold on;

% Tracé des lignes verticales
ylims = get(gca,'YLim');
hold on;
plot([3*H/(0.03) 3*H/(0.03)],ylims,'--','Color','black','LineWidth',3);
text(3*H/(0.03),-0.03, ' \epsilon_{Min,2D}', 'HorizontalAlignment', 'center','FontSize', 28);
hold on;

legend('P1 method','MsFEM dirichlet 3H method','MsFEM Dirichlet 9H method', 'MsFEM periodique 3H method','MsFEM periodique 9H method','FontSize',25)

xlabel("3H/\epsilon",'FontSize',25)
ylabel("Err_{relative}",'FontSize',25)
title(sprintf("Evolution de l'erreur relative H1 en fonction de eps / H=%g",H),'FontSize',25)
set(gca,'FontSize',30)
hold off;
