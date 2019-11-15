function [kab,f] = localk(bnd,x,eNoN,nsd,Ng,J,Wg,ng,fun)


    kab = zeros(eNoN,eNoN);
    f = zeros(eNoN,1);
    xgt = zeros(nsd,1);
%   Loop over one element
    for a = 1:eNoN
        
        if bnd(a)==0
%       Loop over next nodes for kab
        for b = 1:eNoN
            if bnd(b)==0
%           For each node combo, integrate by looping and summing
%           over gauss points
            for g = 1:3
                kab(a,b) = kab(a,b) + Ng(a,g)*Ng(b,g)*J*Wg(g);
            end
            
            elseif bnd(b) ==2
                for g = 1:ng
                    f(a) = f(a) - Ng(a,g)*Ng(b,g)*J*Wg(g);
                end
            end
            
        end
        
%       Now for the f matrix. Loop over gauss points to perform integral
        for g = 1:ng
            
%           To evaluate f(x), we need to get global x at parent gauss pts
            for j = 1:eNoN
                xgt = xgt + Ng(j,g)*x(:,j);
            end
            
%           Now assemble the f matrix by integrating using Gauss Quad
            f(a) = f(a) + Ng(a,g)*fun(xgt(1),xgt(2))*J*Wg(g);
            
%           Reset temp Gauss pts
            xgt(:) = 0;
        end
        
        end
        
    end



end
