function [KAB,F] = globalk(kab,f,nt,IEN)

KAB = zeros(nt,nt);
F = zeros(nt,1);

KAB(IEN,IEN) = kab;
F(IEN) = f;
end