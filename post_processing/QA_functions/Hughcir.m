function[I2,center]=Hughcir(I)

[r c]=size(I);

%defining max &min values for radius
rmax=272;
rmin=271;

PL=zeros(r,c,rmax-rmin+1);
for rad=rmin:rmax
    
    k=1;
    
    % creating a circle for given radius
    for theta=1:360
        x=rad*cosd(theta);
        y=rad*sind(theta);
        X(k)=floor(x+.5);
        Y(k)=floor(y+.5);
        k=k+1;
    end
    
    for i=1:3:r
        for j=1:3:c
            
            if(I(i,j)==1)
                Xi=X+i;
                Yj=Y+j;
                
                index=find((Xi>=rad)&(Yj>=rad)&(Xi<=r-rad)&(Yj<=c-rad));
                
                [rr cc]=size(index);
                P1=logical(zeros(r,c));   
                for l=1:cc
                    P1(Xi(index(l)),Yj(index(l)))=1;
                end
                PL(:,:,rad-rmin+1)=PL(:,:,rad-rmin+1)+P1(:,:);
            end
            
            
        end
    end
end

ma=max(max(max(PL)));

[r c1 p]=find(PL==ma);

pag=floor(c1/c);
ac=c1-(c*pag);

I2=logical(zeros(size(I)));
nrad=pag+rmin;
center(1)=r(1);
center(2)=ac(1);
for theta=1:.2:360
    x=nrad(1)*cosd(theta);
    y=nrad(1)*sind(theta);
    Xn=r(1)+floor(x+.5);
    Yn=ac(1)+floor(y+.5);
    I2(Xn,Yn)=1;
end
