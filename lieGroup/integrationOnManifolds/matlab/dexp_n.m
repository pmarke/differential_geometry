function de = dexp_n(v)

adp = adx(v);
adn = adx(-v);

de = (eye(6)-expm(adn))*inv(adp);

end