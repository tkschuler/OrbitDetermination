
function [E] = anom_ecc(M,e)
    % function [E] = anom_ecc(M,e) 
    % Risoluzione numerica dell'equazione: E-e*sin(E)=M
    % E = anomalia eccentrica [rad]
    % e = eccentricit�
    % M = anomalia media [rad]
    % Si devono dare in input alla funzione due valori scalari,
    % rispettivamente M in rad ed e.
    % Il programma restituisce E [rad] con un errore di 1e-10
    % N.B. Per migliorare l'accuratezza accedere all'editor e modificare il
    % valore della variabile err, che rappresenta l'errore commesso.
    
    format long
%     x = 0;
%     sx = 1;
%     dymax = -1;
%     trovato = true;
%     while (trovato)
%         if (sx<0.2)
%             sx = 0.1 - (x/1000);
%         else
%             sx  = M-x;
%             dx  = M+x;
%             dy  = abs(cos(dx));
%             dy2 = abs(cos(sx));
%         end
%         if (dymax<dy || dymax<dy2)
%             if (dy<dy2)
%                 dymax = dy2;
%             else
%                 dymax = dy;
%                 dy = dymax;
%             end
%         end 
%         f0 = sx-e.*sin(sx)-M;
%         f1 = dx-e.*sin(dx)-M;
%         trovato = ((f0.*f1)>0);
%         x = x + 0.1;
%     end
    E = M;
    k = 1;
    err = 1e-10;
    % stabilito intervallo di ricerca (sx,dx) e max valore derivata dymax;
    while (k>err)
        y = e*sin(E)+M;
        k = abs(abs(E)-abs(y));
        E = y;
    end
    % trovato E con un errore inferiore a 0.1*10^(-4);
    %fprintf(' La soluzione E � stata trovata nell''intervallo [%2.5f,%2.5f]',sx,dx);
    %fprintf(' \n errore inferiore a %1.2e: [rad] E = %3.5f',err,E);
end
