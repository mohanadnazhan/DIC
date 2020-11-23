# DIC
DIC Code / Plate Without hole – Code 1
location = strcat ('D:\Main Project\Test\Aluminum without hole\DSC_00',strcat (num2str(21),'.jpg'));
orgpic= imread(location);
orgpic = rgb2gray (orgpic);
defpic = imread(strcat ('D:\Main Project\Test\Aluminum without hole\DSC_00',strcat (num2str(32),'.jpg')));
defpic= rgb2gray (defpic);
widthpic = size (orgpic,1);
hightpic =size (orgpic,2);
 
 
filter = input ('enter the filter: ')
 
x = filter*2+2:1000:widthpic-filter*2+2;
y = filter*2+2:1000:hightpic-filter*2+2;
checkingvaluex = size (x);
checkingvaluey = size (y);
 
if (checkingvaluex (1,2) < checkingvaluey (1,2))
        zzzz = checkingvaluex(1,2);
else
        zzzz = checkingvaluey(1,2);
end
opopop = input ('Number of square that you want to test: ');
    if ( opopop <= zzzz )
        zzzz=opopop ;
    end
xpoints = zeros(1,zzzz);
u =zeros(1,zzzz);
ypointshifting=zeros (1,zzzz);
v=zeros(1,zzzz);
xp = zeros (1,zzzz);
yp = zeros (1,zzzz);  
 
 
k=1;
k1= 1;
while (k <=zzzz)
t=1;
while (t>0)
% x= input ('x point for the center: ');
% y= input ('FFy point for the center: ');
 
for popoy = 0:15:16
    for popox = 0:15:16
roipic = zeros ((filter*2+1));
z1=1;	
for i=x(k)+popox-filter:x(k)+filter+popox
    z2=1;
    for j=y(k)-filter+popoy:y(k)+filter+popoy
        roipic(z1,z2)=orgpic (i,j);
        z2=z2+1;;
    end 
    z1=z1+1;
end
t=0;
xp (k1)= x(k)+popox;
yp (k1)= y(k)+popoy;
[ypeak xpeak] = dic(defpic,roipic)
yoffSet = ypeak;
xoffSet = xpeak;
figure
imshow(defpic);
imrect(gca, [xp(k1), yp(k1), size(roipic,2), size(roipic,1)]);
shiftingxpoint = x(k)+popox - xpeak-size (roipic,1);
shiftingypoint = y(k)+popoy - ypeak-size (roipic,2);
%%%%%%% 
xpoints(k1) = xpeak;
u(k1) = shiftingxpoint;
ypoints (k1) = ypeak;
v(k1)= shiftingypoint;
k1= k1+1;
end
end
k=k+1;
end
end
 
exi=1/sqrt(3); eta =1/sqrt (3);
ym= 2e11;
por=0.3;
t=1;
% 4 node cartizine Coordinates
%xpoints (1)= 0 ; ypoints (1)= 0;
%xpoints (2) =0.4 ; ypoints (2)= 0;
%xpoints (3)=0.4 ; ypoints (3)=0.3 ;
%xpoints (4)=0 ; ypoints (4)=0.3 ;
% D Matrix
D1= (ym/(1-por^2)) ; D2= (ym*por)/(1-por^2); D3=0; D4=ym/2;
D=[D1 D2 D3; D2 D1 D3 ; D3 D3 D4 ];
%Shape Functions and thier derivates
N1=0.25*(1-exi)*(1-eta) ; N2= 0.25*(1+exi)*(1-eta); N3= 0.25*(1+exi)*(1+eta); N4=0.25*(1-exi)*(1+eta);
dN1dxi= -(1-eta)/4 ; dN2dxi= (1-eta)/4; dN3dxi= (1+eta)/4 ; dN4dxi= -(1+eta)/4;
dN1deta= -(1-exi)/4;dN2deta= -(1+exi)/4; dN3deta= (1+exi)/4; dN4deta = (1-exi)/4;
i1=1;
i2=2;
i3=3;
i4=4;
straindisp= zeros (3,zzzz);
stressdisp = zeros (3,zzzz);
for wk = 1:zzzz
% J Matrix
J11= xpoints (i1)*dN1dxi + xpoints (i2)*dN2dxi+ xpoints (i3)* dN3dxi+xpoints (i4)*dN4dxi;
J12=ypoints (i1)*dN1dxi+ypoints (i2)*dN2dxi+ypoints (i3)*dN3dxi+ypoints (i4)*dN4dxi;
J21=xpoints (i1)*dN1deta+xpoints (i2)*dN2deta+xpoints (i3)*dN3deta+xpoints (i4)*dN4deta;
J22=ypoints (i1)*dN1deta+ypoints (i2)*dN2deta+ypoints (i3)*dN3deta+ypoints (i4)*dN4deta;
J=[ J11 J12; J21 J22];
Jdet= det(J);
% b1 Matrix
B1=[1 0 0 0 ; 0 0 0 1; 0 1 1 0];
B2=[J22/Jdet -J12/Jdet 0 0 ; -J21/Jdet J11/Jdet 0 0 ; 0 0 J22/Jdet -J12/Jdet; 0 0 -J21/Jdet J11/Jdet];
B3=[dN1dxi 0 dN2dxi 0 dN3dxi 0 dN4dxi 0 ; dN1deta 0 dN2deta 0 dN3deta 0 dN4deta 0;
    0 dN1dxi 0 dN2dxi 0 dN3dxi 0 dN4dxi;  0 dN1deta 0 dN2deta 0 dN3deta 0 dN4deta];
B=B1*B2*B3,
% K Matrix
BT=transpose(B);
K = BT*D*B*t*Jdet;
displacement=[u(i1);v(i1);u(i2);v(i2);u(i3);v(i3);u(i4);v(i4)];
Strain=B*displacement;
Stress=D*Strain;
i1=i1+4;
i2=i2+4;
i3=i3+4;
i4=i4+4;
straindisp (:,wk) = Strain;
stressdisp (:,wk) = Stress;
end
 
%%%%%%%%%%%%%%%%%%%%%% location in plates 
locationx =zeros (1,zzzz);
locationy =zeros (1,zzzz);
locationxpeak =zeros (1,zzzz);
locationypeak =zeros (1,zzzz);
locationxshift =zeros (1,zzzz);
locationxshift =zeros (1,zzzz);
orgimg = size (orgpic);
defimg = size (defpic);
twidth  = input ('enter the width of plate'); 
thight = input ('enter the hight of palte');
for i=1:zzzz
locationx (i) = xp (i) / orgimg (1) * twidth ;
locationy (i) = yp (i) / orgimg (2) * thight ;
locationxpeak (i) = double (xpoints (i) )/ double (defpic (1)) * twidth ;
locationypeak (i)= double (ypoints (i)) /double( defpic (2) )* thight ;
locationxshift (i) =double(u (i) )/ double (defpic(1) )* twidth *10^-3 ;
locationxshift (i)=double(v (i) )/ double(defpic (2) )* thight ;
end
 
%%%% results 
locationx 
locationy 
locationxpeak 
locationypeak 
locationxshift
locationxshift
 
 
%%%%% eq
eqx =fix (sum (xpoints)/4)
eqy = fix(sum (ypoints)/4)
 
% plot (straindisp,stressdisp)
dsiplay (straindisp ,stressdisp,displacement)
