function val = condition_initiale(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% condition_initiale :
% Evaluation de la condition initiale.
%
% SYNOPSIS val = condition_initiale(x,y)
%          
% INPUT * x : coordonnee du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%val = (10+15*(abs(x-0.4)<0.25)).*sin(pi*x);

%val = 1.*sin(pi*x);

%val = (1.*sin(3*pi*x) + sin(3*sqrt(2)*pi*x) + x +2).*x.*(1-x); %Condition
%initiale plutot simple

%val = (1.*sin(12*pi*x) + sin(12*sqrt(2)*pi*x) + x +2).*x.*(1-x);%Condition

val = 3*(1.*sin(7*pi*x) + sin(7*sqrt(2)*pi*x) + x +2).*x.*(1-x).*(1-x).*(1-x);%Condition
%initiale plus compliquée
%val = 10.*sin(10*pi*x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
