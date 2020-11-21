u1 = rand(1000,1);
u2 = rand(1000,1);



a = sqrt(-2*log(u1)).*cos(2*pi*u2);

% for ii = 1:length(a)
% 
%     b = rand(1,1);
%     if (b <= 0.5)
%         a(ii) = a(ii) *(-1);
%     end
%     
% end




figure(1),clf;
histogram(a,'Normalization','probability')



function p = myNormal(x)

1/sqrt(2*pi)*exp(-1/2*x^2);

end