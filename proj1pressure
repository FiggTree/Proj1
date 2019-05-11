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
    a=min(a);
    b=min(b);
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
hold('on')
plot(H2Omolamt(v(1,2))./10,d,'ro')
xlabel('Mole fraction H2O put in')
ylabel('Mole fraction of H2 out')
legend('H2vsH2O','maximum')
title('Mole percent H20 put in versus Mole Fraction H2 produced')


%figure
%plot(H2Omolamt./10, extent(:,1),'r')
%hold on
%plot(H2Omolamt./10, extent(:,2),'b')
%plot(H2Omolamt./10, (3.*extent(:,1)+extent(:,2)),'k:')
%plot(H2Omolamt(f)./10,(3.*extent(f,1)+extent(f,2)), 'kx')


PercentH2O = (H2Omolamt(v(1,2))./10);
PercentCH4 = 1-PercentH2O;
FractionH2 = v(1,1);
extentr1= extent(f,1);
extentr2= extent(f,2);
l=['Using the initial condidtions with the fraction of H2O being ', num2str(PercentH2O), ' and the fraction of CH4 being ', num2str(PercentCH4), ' gives the maximum fraction of H2 out: ', num2str(FractionH2),'. The extents of reaction at this time are e1=', num2str(extentr1),' and e2=', num2str(extentr2),'.'];
disp(l) 
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
    k1=((((Yco(e1,e2)).*((Yh2(e1,e2)).^3))./(((Ych4(e1,e2,(10-molh2o))).*(Yh2o(e1,e2,molh2o)))))-(1.93.*10.^(-4)));
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
