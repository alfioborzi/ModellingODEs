function [ v ] = argminH(y,u,p,A,B,epsilon,OCP )

 


    num_Int=100;
    Disc_KU=linspace(OCP.umin,OCP.umax,num_Int);
    dist=inf;

    while(dist>10^-14)
    [~,pos]=min((OCP.nu/2)*Disc_KU.^2+OCP.beta*abs(Disc_KU) + (p')*A*y+(p')*B*y.*Disc_KU+epsilon*((Disc_KU-u).^2));
    v=Disc_KU(pos);

    if (pos==1)
    Disc_KU=linspace(Disc_KU(pos),Disc_KU(pos+1),num_Int);
    dist=(Disc_KU(pos+1)-Disc_KU(pos));
    elseif  (pos==max(size(Disc_KU)))
    Disc_KU=linspace(Disc_KU(pos-1),Disc_KU(pos),num_Int);
    dist=(Disc_KU(pos)-Disc_KU(pos-1));
    else
    Disc_KU=linspace(Disc_KU(pos-1),Disc_KU(pos+1),num_Int);
    dist=Disc_KU(pos+1)-Disc_KU(pos-1);
    end
    end

end


