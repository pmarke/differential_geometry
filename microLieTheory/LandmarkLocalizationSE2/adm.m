function ad = adm(tau)

ad = [ssm(tau(3)), -ssm(1)*tau(1:2); zeros(1,3)];

end