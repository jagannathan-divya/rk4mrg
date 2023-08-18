function [expLh, M1, M2] = expMat(L,dt)
   [V,D] = eig(dt*L);
   expLh = (V*expm(D))/(V);
 
   M1 = L\(expLh-eye(size(L)));
   L2 = L*L;
   M2 = L2\(expLh-(eye(size(L))+dt*L));
end