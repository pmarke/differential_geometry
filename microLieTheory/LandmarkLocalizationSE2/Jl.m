function J = Jl(tau)


v1 = tau(1);
v2 = tau(2);
th = tau(3);


ad = adm(tau);

if (th < 1e-4)
    
    J = zeros(3,3);
    
    for ii = 0:10
        J = J + (ad)^ii/factorial(ii+1);
    end

else

J = [sin(th)/th,     (cos(th)-1)/th,   (th*v1+v2-v2*cos(th)-v1*sin(th))/th^2;...
      (1-cos(th))/th, sin(th)/th,       (-v1+th*v2+v1*cos(th)-v2*sin(th))/th^2;...
      0 0 1];
end



end
