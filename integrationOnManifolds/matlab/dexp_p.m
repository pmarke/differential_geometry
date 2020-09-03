function de = dexp_p(v)

ad = adx(v);

de = inv(ad)*(eye(6)-expm(ad))
de = (expm(ad)-eye(6))*inv(ad)

end