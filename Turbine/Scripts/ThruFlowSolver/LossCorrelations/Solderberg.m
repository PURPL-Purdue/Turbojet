function xi = Solderberg(alpha_1, alpha_2, AR, t2l, Re)
    A = [5.35647E-08 -2.96536E-06 0.000358208 0.042072023];
    B = [-0.0948064 0.019066743];

    eps = rad2deg(alpha_1)+rad2deg(alpha_2);
    xi_st = dot(A, [eps^3 eps^2 eps 1]) + dot(B, [t2l 1]);
    xi_pr = (1+xi_st)*(0.975+0.075/AR)-1;
    xi = xi_pr * (10^5/Re)^0.25;
end