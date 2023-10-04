
function dfx = dropdfdxXY(xj,yj,tj,p) 

dfx = (-sin(tj))/(xj*yj+p*xj-sin(tj))^2; 

   