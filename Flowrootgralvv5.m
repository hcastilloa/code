clear
clc
%INTERVAL Ci AND Cij:
Liminfa=0.0001;
Liminf=0;
LimSup=8;
d=0.0125;
e=1./d;
k=(LimSup./d).^2;  %Number of rows%
%MATRIX Cij%
aa1=1.90;
aa2=100000.00;
s=linspace(Liminfa,LimSup,LimSup./d);
H=s;
for i=1:(LimSup./d)-1
    C12=[s,H];
    H=C12;
end
count16=1;
for i=1:k
       if H(1,i)<LimSup
       C21(1,i)=H(1,count16);
       else
           C21(1,i)=H(1,count16);
           count16=count16+1;
          end
end
for i=1:k
    C1(i,1)=aa1;
end
for i=1:k
    C2(i,1)=aa2;
end
format short
Coefac=[C1,C2,H',C21'];
%_______________________________________________________________________________________________________________%
%COEFFICIENT GAMMA OF GENERAL POLYNOMIAL G(Ci,Cij):
Ga4=-(Coefac(:,1).^2).*Coefac(:,2);
Ga3=((2.*(Coefac(:,1)-1)).*Coefac(:,1).*Coefac(:,2));
Ga2=-Coefac(:,2).*(Coefac(:,1)-1).^2-Coefac(:,1).*(Coefac(:,2)-1).*Coefac(:,3)+2.*Coefac(:,2).*Coefac(:,1);
Ga1=-2.*Coefac(:,2).*(Coefac(:,1)-1)-Coefac(:,4).*Coefac(:,3).^2+(Coefac(:,1)-1).*(Coefac(:,2)-1).*Coefac(:,3);
Ga0=(Coefac(:,2)+Coefac(:,3)).*(Coefac(:,3)-1);
Mgamma=[Ga4,Ga3,Ga2,Ga1,Ga0];
%___________________________________________________________________________________________________________________%
%MATRIZ DE CEROS DEL POLONOMIO GENERAL EN COMPETICION POR INT:
for j= 1:k
Pn=Mgamma(j,:);
RootPn=roots(Pn);
S(j,:)=RootPn;
end
format long
S;

% % %________________________________________________________________________________________________________________%
%POSITIVE ROOTS OF POLINOMIAL
[z,o]=size(S);
count=0;
for i=1:k
    count25=0;
        for n=1:o
           if imag(S(i,n))==0
             if real(S(i,n))>0
               count25=count25+1;
               count=count+1;
               CS(i,1)=count25; 
             LPRTri(i,n)=S(i,n); %Puntos U1 para Tri-estabilidad            
             end
           
           end
        end
           end
format short
%**********************************************Tri-stability*****************************************************%
     LPRTricount=[LPRTri,Coefac,CS];
     [wq,en]=size(LPRTricount);
     count22=0;
     count23=0;
     count24=0;
     %one Root 
     for i=1:wq
     if LPRTricount(i,en)==1
         count22=count22+1;
         LPRTricount1(count22,:)=LPRTricount(i,:);       
     end    
     %Two Roots
     if LPRTricount(i,en)==2
         count23=count23+1;
         LPRTricount2(count23,:)=LPRTricount(i,:);         
     end
     %Three Roots
     if  LPRTricount(i,en)==3
         count24=count24+1;
         LPRTricount3(count24,:)=LPRTricount(i,:);  
       end
     end 
     
     %COLUMNS ROOTS FOURTH DEGREE (WE MUST IMPROVE IT).
     if aa1.*aa2==0
LTRCcol1=[LPRTricount1(:,1),LPRTricount1(:,3:7);LPRTricount1(:,2),LPRTricount1(:,3:7)];     
LTRCcol2=[LPRTricount2(:,1),LPRTricount2(:,3:7);LPRTricount2(:,2),LPRTricount2(:,3:7)];
LTRCcol=[LTRCcol1;LTRCcol2];
     else
 LTRCcol1=[LPRTricount1(:,1),LPRTricount1(:,5:9);LPRTricount1(:,2),LPRTricount1(:,5:9);LPRTricount1(:,3),LPRTricount1(:,5:9);LPRTricount1(:,4),LPRTricount1(:,5:9)];
 LTRCcol2=[LPRTricount2(:,1),LPRTricount2(:,5:9);LPRTricount2(:,2),LPRTricount2(:,5:9);LPRTricount2(:,3),LPRTricount2(:,5:9);LPRTricount2(:,4),LPRTricount2(:,5:9)];
 LTRCcol3=[LPRTricount3(:,1),LPRTricount3(:,5:9);LPRTricount3(:,2),LPRTricount3(:,5:9);LPRTricount3(:,3),LPRTricount3(:,5:9);LPRTricount3(:,4),LPRTricount3(:,5:9)];
 LTRCcol=[LTRCcol1;LTRCcol2;LTRCcol3];
     end
 [cw,mp]=size(LTRCcol);
count28=0;
for i=1:cw
    if LTRCcol(i,1)>0
        count28=count28+1;
        LTRCC(count28,:)=LTRCcol(i,:);
    end
end

%******************************************************************************************************************%     
%Function Nullclines (U1,U2)%
%________________________________________________________________________________________________________________%
[ap,pe]=size(LTRCC);
for i=1:ap
    POR1(i,1)=LTRCC(i,1);
    POR2(i,1)=((1-POR1(i,1)).*(1+LTRCC(i,2).*POR1(i,1)))./(LTRCC(i,4)); 
end
LTRCU=[POR2,LTRCC];
count18=0;
for i=1:ap
    if LTRCU(i,1)>0
            count18=count18+1;
        PORRC(count18,:)=LTRCU(i,:);
        end
end
Uscoef=PORRC;
[er,ty]=size(Uscoef);
%-------Jacobian Matrix-----------%%%%%%
syms u1 u2
for i=1:er  
cc1=Uscoef(i,3);
cc2=Uscoef(i,4);
cc12=Uscoef(i,5);
cc21=Uscoef(i,6);
F=[u1-u1.^2-((cc12.*u1.*u2)./(1+cc1.*u1)),u2-u2.^2-((cc21.*u1.*u2)./(1+cc2.*u2))];
v=[u1 u2];
J=jacobian(F,v);
PP=[Uscoef(i,2) Uscoef(i,1)];
JN=double(subs(J,v,PP));
Lanc(i,:)=eig(JN)';
end
Lan=[Lanc,Uscoef];
% %FORCE AND ABILITY TO LIMIT SPECIES TO STABILITY%
% %________________________________________________________________________________________________________________%
[qw,as]=size(Lan);
count3=0;
for i=1:qw
        if Lan(i,1)<0
            if Lan(i,2)<0
                count3=count3+1;
                Lanneg(count3,:)=Lan(i,:);
            end
        end     
end
[m,n]=size(Lanneg);
%graphics Cij%
%________________________________________________________________________________________________________________%
count4=0;
for i=1:m
if Lanneg(i,7)<1
    if Lanneg(i,8)<1
        count4=count4+1;
        coe1(count4,1)=Lanneg(i,7);
        coe2(count4,:)=Lanneg(i,[8 9]);
        coe=[coe1,coe2];
    end
end
end
count11=0;
for i=1:m
if Lanneg(i,7)>1
    if Lanneg(i,8)>1
        
        count11=count11+1;
        exclu1(count11,1)=Lanneg(i,7);
        exclu2(count11,:)=Lanneg(i,[8 9]);
        exclu=[exclu1,exclu2];
    end
end
end


count12=0;
for i=1:m
if Lanneg(i,7)>=1
    if Lanneg(i,8)<=1
        count12=count12+1;
        sp21(count12,1)=Lanneg(i,7);
        sp22(count12,1)=Lanneg(i,8);
        sp2=[sp21,sp22];
    end
end
end


count13=0;
for i=1:m
if Lanneg(i,7)<1
    if Lanneg(i,8)>1
        count13=count13+1;
        sp11(count13,1)=Lanneg(i,7);
        sp12(count13,1)=Lanneg(i,8);
          sp1=[sp11,sp12];
    end
end
end
figure
if count4>0
scatter(coe(:,1),coe(:,2),'.','g')
end
hold on
if count11>0
scatter(exclu(:,1),exclu(:,2),'+','k')
end
if count13>0
scatter(sp1(:,1),sp1(:,2),'*','r')
end
if count12>0
scatter(sp2(:,1),sp2(:,2),'x','b')
end
% ylim([0 2])
% xlim([0 2])
xlabel('c_{12}'),ylabel('c_{21}')
hold off
grid on
grid minor
