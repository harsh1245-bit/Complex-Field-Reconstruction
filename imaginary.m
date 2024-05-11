clear all;
close all;
clc;
rows = 400;
cols = 400;
squareSize = 4;
% Initialize the checkerboard matrix
board = zeros(rows, cols);
% Iterate through rows and columns to set values
for i = 1:rows
    for j = 1:cols
        % Calculate the square index for the current cell
        squareRow = ceil(i / squareSize);
        squareCol = ceil(j / squareSize);
        
        % Check if the sum of squareRow and squareCol is odd
        if mod(squareRow + squareCol, 2) == 1
            % If odd, set the value to 1
            board(i, j) = 1;
        end
    end
end
% Parameters
period = 4;           % Period of the grating in pixels
dimensions = 400;      % Dimensions of the square grating in pixels
% Create meshgrid
[x, y] = meshgrid(1:dimensions, 1:dimensions);
% Calculate grating using a 45-degree rotation
grating = mod(floor(x * cos(pi/4)+y*sin(pi/4) ), period) < period / 2;
% Display the grating
figure(1);
imagesc(grating);
colormap(gray);
title('45-degree Rotated Binary Grating');
axis equal;
axis tight;
lambda=632.8*10^(-9);
k=2*pi/lambda;
w=0.001;
l=1;
L=0.01;
n=399;
dx=L/n;
x=(-L/2:dx:L/2);
[X,Y] = meshgrid(x,x);
r=sqrt(X.^2+Y.^2);
phase=atan2(Y,X);
 f=exp(1i*phase*l);
 figure(2);imagesc(angle(f));
 
%imwrite(f,'C:\Users\HP\Desktop\Holograms\l=0.jpg','jpg') 
E=exp(-r.^2/w^2).*(X+sign(l).*1i*Y).^abs(l); %.*((X+1i*Y).^l);%.*exp(i*k*(r^2)/2*R);
figure(3); imagesc(abs(E));
title('vortex');
comb = grating.*E;
figure(4); imagesc(abs(comb));
for_0 = comb;
for_pi = comb;
for i = 1:400
   for j = 1:400
      if board(i,j)==1
          for_0(i,j) = comb(i,j).*exp(1i * pi/2);
	for_pi(i,j) = comb(i,j).*exp(1i * 3*pi/2);
	end
   end
end
figure(5); imagesc(angle(for_0));
figure(6); imagesc(angle(for_pi));
%%% sinusoidal patterns 
a=1;
b=1;
lx1=2.12;
ly1=lx1;
dx1=5.3e-3;
dy1=dx1;
[m1,n1]=meshgrid(-0.5*lx1:dx1:0.5*lx1,+0.5*ly1:-dy1:-0.5*ly1);
for x=1:400
         for y=1:400
             T(x,y)=(2*pi*0.02*(x-201));
             R(x,y)=(2*pi*0.02*(y-201));
     
        end
end
count=0;
for l=1:10:400
        for p=1:10:400
            %fprintf('grnd: %d\n', vv)
            fprintf('Iteration: %d\n', l)
            fprintf('sub: %d\n', p)
            tic
            p_0=a+b.*cos((m1(l,p)*T)+(n1(l,p)*R)+0);
            p_pi= a+b.*cos((m1(l,p)*T)+(n1(l,p)*R)+pi);
            p_pibi2=a+b.*cos((m1(l,p)*T)+(n1(l,p)*R)+pi/2);
            p_3pibi2=a+b.*cos((m1(l,p)*T)+(n1(l,p)*R)+3*pi/2);
            
         
              
                 K_0_0 = (for_0.* p_0);
		 K_0_pi = (for_pi.*p_0);
                 Sum_0_0 = sum(K_0_0(:));
		 Sum_0_pi = sum(K_0_pi(:));
                 I_0_0=(Sum_0_0.*conj(Sum_0_0)); % capturing intensity in CCD 
                 I_0_pi=(Sum_0_pi.*conj(Sum_0_pi)); % capturing intensity in CCD                 
                 I_0 = I_0_0 - I_0_pi;
                  
                 K_pi_0 = (for_0.* p_pi);
		 K_pi_pi = (for_pi.*p_pi);
                 Sum_pi_0 = sum(K_pi_0(:));
		 Sum_pi_pi = sum(K_pi_pi(:));
                 I_pi_0=(Sum_pi_0.*conj(Sum_pi_0)); % capturing intensity in CCD 
                 I_pi_pi=(Sum_pi_pi.*conj(Sum_pi_pi)); % capturing intensity in CCD                 
                 I_pi = I_pi_0 - I_pi_pi;
                 
                 K_pibi2_0 = (for_0.* p_pibi2);
		 K_pibi2_pi = (for_pi.*p_pibi2);
                 Sum_pibi2_0 = sum(K_pibi2_0(:));
		 Sum_pibi2_pi = sum(K_pibi2_pi(:));
                 I_pibi2_0=(Sum_pibi2_0.*conj(Sum_pibi2_0)); % capturing intensity in CCD 
                 I_pibi2_pi=(Sum_pibi2_pi.*conj(Sum_pibi2_pi)); % capturing intensity in CCD                 
                 I_pibi2 = I_pibi2_0 - I_pibi2_pi;
                 
                 K_3pibi2_0 = (for_0.* p_3pibi2);
		 K_3pibi2_pi = (for_pi.*p_3pibi2);
                 Sum_3pibi2_0 = sum(K_3pibi2_0(:));
		 Sum_3pibi2_pi = sum(K_3pibi2_pi(:));
                 I_3pibi2_0=(Sum_3pibi2_0.*conj(Sum_3pibi2_0)); % capturing intensity in CCD 
                 I_3pibi2_pi=(Sum_3pibi2_pi.*conj(Sum_3pibi2_pi)); % capturing intensity in CCD                 
                 I_3pibi2 = I_3pibi2_0 - I_3pibi2_pi;
                 
            
           
            II_0 =0;
            II_pi=0;
            II_pibi2=0;
            II_3pibi2=0;
            
                
               II_0 =  I_0;
                 
               II_pi =  I_pi;
                
               II_pibi2 = I_pibi2;
                
               II_3pibi2 =  I_3pibi2;
                
            
            sin_II_0(l,p)= II_0; % single pixel measurements
     
            Sin_II_pi(l,p)= II_pi;
            
            Sin_II_pibi2(l,p)= II_pibi2;
            
            Sin_II_3pibi2(l,p)= II_3pibi2;
                     
           toc 
            
        end
  end
  
%%%%%%%%
toc
C_xy = (sin_II_0-Sin_II_pi)+1i*(Sin_II_pibi2-Sin_II_3pibi2);
 figure;
 imagesc(abs(C_xy));
 figure;
 tk=fftshift(fft2(C_xy));
 figure;
 imagesc(abs(tk));
 figure;
 imagesc(abs(C_xy));
 figure;
 imagesc(angle(C_xy));
sparse_img=abs(C_xy);
[x, y] = meshgrid(1:391, 1:391);
% define the new grid to interpolate the image on
[Xi, Yi] = meshgrid(1:0.25:391, 1:0.25:391);
% interpolate the image using griddata
interp_img = griddata(x(sparse_img~=0), y(sparse_img~=0), sparse_img(sparse_img~=0), Xi, Yi);
 
figure;
imagesc(interp_img);
