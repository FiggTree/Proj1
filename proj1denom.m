clear
syms e1 e2 positive
n=1;
for i=0:.1:10
    [a b]=vpasolve(K1(e1,e2,i),K2(e1,e2,i));
    if length(a)==0 || length(b)==0
        a=0;
        b=0;
    end
    if i==0 || i==10
        a=0;
        b=0;
    end
    ae1=round(double(a),4);
    be2=round(double(b),4);
    extent(n,1:2)=[ae1 be2];
    CH4molamt(n,1)=10-i;
    H2Omolamt(n,1)=i;
    H2molfrac(n,1)=Yh2(ae1,be2);
    COmolfrac(n,1)=Yco(ae1,be2);
    CO2molfrac(n,1)=Yco2(ae1,be2);
    CH4molfrac(n,1)=Ych4(ae1,be2,10-i);
    H2Omolfrac(n,1)=Yh2o(ae1,be2,i);
    H2molamt(n,1)=(3.*ae1 + be2);
    n=n+1;
end
[d f]=max(H2molfrac);
v=[d f];

figure
plot(H2Omolamt./10,H2molfrac,'k')
hold on
grid on
ylim([0 0.10]);
plot(H2Omolamt(v(1,2))./10,d,'ro','MarkerSize',5)
txt='(0.67, 0.0842)';
text(.52, (H2molfrac(68) +.003), txt)
xlabel('Initial Mole Fraction of H_2O')
ylabel('Mole Fraction of H_2 produced')
legend('H_2 vs H_2O','maximum')
title('Final Mole Fraction H_2 produced versus Initial Mole Fraction H_2O')

figure
plot(H2Omolamt,H2molamt./9,'r',H2Omolamt,H2molfrac,'b')
hold on
ylim([0 0.18]);
xlabel('Mole Fraction of H_2O')
ylabel('H_2 Amount and Fraction')
legend('H_2 Mole Amount (scaled)','H_2 Mole Fraction')
title('Comparison of H_2 Amount and Fraction versus Initial H_2O')

figure
o1=plot(H2Omolamt,extent(:,1));
c1=o1.Color;
o1.Color=[0.6350 0.0780 0.1840];
hold on
o2=plot(H2Omolamt,extent(:,2));
c2=o2.Color;
o2.Color=[0.4940 0.1840 0.5560];
plot(H2Omolamt,2.*((3.*extent(:,1)+extent(:,2))./(10+2.*extent(:,1))),'k--')
ylim([0 .3])
xlabel('Mole Fraction of H_2O')
ylabel('Extents of Reactions and H_2 Mole Fraction')
legend('Extent of Reaction 1', 'Extent of Reaction 2', 'Resulting H_2 Mole Fraction (scaled)')
title('Extents of Reactions and Resulting H_2 Mole Fraction versus Inital H_2O')

PercentH2O = (H2Omolamt(v(1,2))./10);
PercentCH4 = 1-PercentH2O;
FractionH2 = v(1,1);
extentr1= extent(f,1);
extentr2= extent(f,2);
l=['Using the initial condidtions with the fraction of H2O being ', num2str(PercentH2O),' and the fraction of CH4 being ', num2str(PercentCH4), ];
ll=['produces the maximum mole fraction of H2: ', num2str(FractionH2),'.'];
lll=[ 'The extents of reaction at this time are e1=', num2str(extentr1),' and e2=', num2str(extentr2),'.'];
disp(l) 
disp(ll)
disp(lll)
sz=[1,5];
vartypes={'double','double','double','double','double'};
varnames={'CH4','H2O','H2','CO','CO2'};
Equlibrium_Gas_Concentrations=table('Size',sz,'VariableTypes',vartypes,'VariableNames',varnames);
Equlibrium_Gas_Concentrations(1,:)={CH4molfrac(f,1),H2Omolfrac(f,1),H2molfrac(f,1),COmolfrac(f,1),CO2molfrac(f,1)}


function k1 = K1(e1,e2,molh2o)
if molh2o>10 || molh2o<0
    disp('ERROR: mols of h2o exceeds basis range of 0-10')
    return
else
    k1=((((Yco(e1,e2)).*((Yh2(e1,e2)).^3).*((10+2.*e1).^2)./((Ych4(e1,e2,(10-molh2o))).*(Yh2o(e1,e2,molh2o)))))-(1.93.*10.^(-4)));
end
end
function k2 = K2(e1,e2,molh2o)
if molh2o>10 || molh2o<0
    disp('ERROR: mols of h2o exceeds basis range of 0-10')
    return
else
    k2=(((Yco2(e1,e2)).*(Yh2(e1,e2)))./((Yco(e1,e2)).*(Yh2o(e1,e2,molh2o))))-5.528;
end
end
function yh2o = Yh2o(e1,e2,molh2o)
yh2o=((molh2o-e1-e2)./(10+2.*e1));
end
function ych4 = Ych4(e1,e2,molch4)
ych4=((molch4-e1)./(10+2.*e1));
end
function yh2 = Yh2(e1,e2)
yh2=((3.*e1+e2)./(10+2.*e1));
end
function yco = Yco(e1,e2)
yco=((e1-e2)./(10+2.*e1));
end
function yco2 = Yco2(e1,e2)
yco2=(e2./(10+2.*e1));
end
