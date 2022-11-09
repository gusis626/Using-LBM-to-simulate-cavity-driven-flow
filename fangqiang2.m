clear;
clc;

Lx=150; 
Ly=150; 

dL=1; 
dT=1; 
NX=round(Lx/dL);
NY=round(Ly/dL);

% D2Q9
global w;
w =[ 4/9,  1/9,  1/9,  1/9,  1/9, 1/36, 1/36, 1/36, 1/36];
global c;
c =[0, 0; 1, 0; 0, 1;-1, 0; 0,-1; 1, 1;-1, 1;-1,-1; 1,-1];

[x,y]=meshgrid(1:NX+1,1:NY+1);
density = zeros(NY+1,NX+1);
u       = zeros(NY+1,NX+1,2);
f       = zeros(NY+1,NX+1,9);
f_eq    = zeros(NY+1,NX+1,9);

%initialize
c_s =sqrt(1/3)*dL/dT;
density_0=1;  
nu   =0.01;
Re=5000;
U=nu*Re/Lx;
tau=nu/(c_s^2)+0.5*dT;

density(:,:)=density_0;
u(NY+1,:,1)=U;


objvideo=VideoWriter('example9.avi');
objvideo.FrameRate=round(20);
open(objvideo)
H=525*3;
W=700*3;

MaxIteration=70000;
for ite=1:MaxIteration
    
    % collision
    for k=1:9
        f_eq=w(k)*density(:,:).*(1+(u(:,:,1)*c(k,1)+u(:,:,2)*c(k,2))/(c_s^2)+((u(:,:,1)*c(k,1)+u(:,:,2)*c(k,2)).^2)/(2*c_s^4)-(u(:,:,1).*u(:,:,1)+u(:,:,2).*u(:,:,2))/(2*c_s^2));
        f(:,:,k)=f(:,:,k)*(1-dT/tau)+f_eq*dT/tau;
    end
    
    % streaming
    f(:,2:NX+1,2) = f(:,  1:NX,2);
    f(:,  1:NX,4) = f(:,2:NX+1,4);
    f(2:NY+1,:,3) = f(1:NY  ,:,3);
    f(1:NY  ,:,5) = f(2:NY+1,:,5);
    f(2:NY+1,2:NX+1,6) = f(  1:NY,  1:NX,6); 
    f(2:NY+1,  1:NX,7) = f(  1:NY,2:NX+1,7); 
    f(1:NY  ,1:NX  ,8) = f(2:NY+1,2:NX+1,8); 
    f(1:NY  ,2:NX+1,9) = f(2:NY+1,  1:NX,9); 
    
    % boundary condition
    f(2:NY+1,1,2)=f(2:NY+1,1,4);
    f(2:NY+1,1,6)=f(2:NY+1,1,8);
    f(2:NY+1,1,9)=f(2:NY+1,1,7);
    
    f(2:NY+1,NX+1,4)=f(2:NY+1,NX+1,2);
    f(2:NY+1,NX+1,7)=f(2:NY+1,NX+1,9);
    f(2:NY+1,NX+1,8)=f(2:NY+1,NX+1,6);
    
    f(1,:,3)=f(1,:,5);
    f(1,:,6)=f(1,:,8);
    f(1,:,7)=f(1,:,9);
    
    density_temp=f(NY+1,2:NX,1)+f(NY+1,2:NX,2)+f(NY+1,2:NX,4)+2*(f(NY+1,2:NX,3)+f(NY+1,2:NX,6)+f(NY+1,2:NX,7));
    f(NY+1,2:NX,5)=f(NY+1,2:NX,3);
    f(NY+1,2:NX,9)=f(NY+1,2:NX,7)+ density_temp * U /6;
    f(NY+1,2:NX,8)=f(NY+1,2:NX,6)- density_temp * U /6;
    
    % density
    density=sum(f,3);
    
    %velocity
    u=zeros(NY+1,NX+1,2);
    for k=1:9
        u(:,:,1)=u(:,:,1)+f(:,:,k)*c(k,1);
        u(:,:,2)=u(:,:,2)+f(:,:,k)*c(k,2);
    end
    u=u./density;
    
    u(NY+1,2:NX,1)=U;
    u(NY+1,2:NX,2)=0;
    
    if mod(ite,100)==0
    u_norm=sqrt(u(:,:,1).^2+u(:,:,2).^2);
    str=zeros(150,150);
for i=1:150
for j=2:150
str(i,j)=str(i,j-1)+0.25*(density(i,j)+density(i,j-1))*(u(i,j,2)+u(i,j-1,2));
end
end
startx=0:15:150;
starty=0:15:150;
[startx1,starty1]=meshgrid(startx,starty);
x=1:151;
y=1:151;

    frame=getframe(gcf);
    frame.cdata=imresize(frame.cdata,[H W]);
subplot(2,2,1)
    pcolor(x,y,u_norm);
    shading interp;
    hcb=colorbar;
    title(ite);
    axis equal;
    view([0,0,1]);
    subplot(2,2,2)
    contour(str,100)
    axis equal
    subplot(2,2,3)
yy=plot(0,0)
delete(yy)
  yy= streamline(x,y,u(:,:,1),u(:,:,2),startx1,starty1)
axis equal
    title('RE=5000')
    axis equal;

    writeVideo(objvideo,frame)
end
end

close(objvideo);