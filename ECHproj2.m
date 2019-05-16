clear
syms A B V Cwd Cs Cw E Gwd Gs Gw Gg Pwd Ps Pw Pg
eq1= B+A==1000;
eq2= .65.*A==Cwd;
eq3= (1./3).*.3.*A==Cw;
eq4= .05.*A==Cs;
eq5= (2./3).*.3.*A==V;
eq6= Cwd==Gwd;
eq7= Cw==Gw;
eq8= Cs==Gs;
eq9= E==Gg;
eq10= E==.5.*Cw;
eq11= Gwd+.65.*B==Pwd;
eq12= Gw+.3.*B==Pw;
eq13= Gs+.05.*B==Ps;
eq14= Gg==Pg;
eq15= Pg./(Pwd+Pw+Ps+Pg)==.03;

a=solve(eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13,eq14,eq15,A,B,V,Cwd,Cs,Cw,E,Gwd,Gs,Gw,Gg,Pwd,Ps,Pw,Pg)
