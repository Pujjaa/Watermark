clc;
clear;
close all;
cover=double(imread('D:\project 17\171.jpeg'));
cover=imresize(cover,[512,512]);
cover_red=cover(:,:,1);
cover_green=cover(:,:,2);
cover_blue=cover(:,:,3);
watermarked=zeros(512,512,3);
watermark=imread('D:\project 17\172.jpeg');
watermark=imresize(watermark,[64,64]);
watermark_red=watermark(:,:,1);
watermark_green=watermark(:,:,2);
watermark_blue=watermark(:,:,3);
H_N=[1 1 1 1; 1 -1 1 -1;1 1 -1 -1;1 -1 -1 1];
%%%%%%%%%%%%%%%%%%%%%%%%% ARNOLD TRANSFORM %%%%%%%%%%%%%%%%%%%
watermark_scramble=zeros(64,64);
watermark_extracted=zeros(64,64);
watermark_rescramble=zeros(64,64);

watermark_scramble_red=zeros(64,64);
watermark_extracted_red=zeros(64,64);
watermark_rescramble_red=zeros(64,64);

watermark_scramble_green=zeros(64,64);
watermark_extracted_green=zeros(64,64);
watermark_rescramble_green=zeros(64,64);

watermark_scramble_blue=zeros(64,64);
watermark_extracted_blue=zeros(64,64);
watermark_rescramble_blue=zeros(64,64);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    SRAMBLING    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=0:63
    for j=0:63
        coordinator=[1,1;1,2]*[i;j];
        m=rem(coordinator(1),64);
        n=rem(coordinator(2),64);
        watermark_scramble_red(m+1,n+1)=watermark_red(i+1,j+1);
    end
end

for i=0:63
    for j=0:63
        coordinator=[1,1;1,2]*[i;j];
        m=rem(coordinator(1),64);
        n=rem(coordinator(2),64);
        watermark_scramble_green(m+1,n+1)=watermark_green(i+1,j+1);
    end
end

for i=0:63
    for j=0:63
        coordinator=[1,1;1,2]*[i;j];
        m=rem(coordinator(1),64);
        n=rem(coordinator(2),64);
        watermark_scramble_blue(m+1,n+1)=watermark_blue(i+1,j+1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Decimal to binary %%%%%%%%%%%%%%%%%%%%%%%%%
watermark_scramble_red=reshape(watermark_scramble_red,[64*64,1]);
watermark_scramble_green=reshape(watermark_scramble_green,[64*64,1]);
watermark_scramble_blue=reshape(watermark_scramble_blue,[64*64,1]);
watermark_scramble_red=de2bi(watermark_scramble_red,8);
watermark_scramble_green=de2bi(watermark_scramble_green,8);
watermark_scramble_blue=de2bi(watermark_scramble_blue,8);
watermark_scramble_red=reshape(watermark_scramble_red,[64*64*8,1]);
watermark_scramble_green=reshape(watermark_scramble_green,[64*64*8,1]);
watermark_scramble_blue=reshape(watermark_scramble_blue,[64*64*8,1]);
watermark_scramble_red_extract=watermark_scramble_red;
watermark_scramble_green_extract=watermark_scramble_green;
watermark_scramble_blue_extract=watermark_scramble_blue;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EMBEDDING %%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:128
    for j =1:128
        redblock=cover_red(4*(i-1)+1:4*i,4*(j-1)+1:4*j);
        greenblock=cover_green(4*(i-1)+1:4*i,4*(j-1)+1:4*j);
        blueblock=cover_blue(4*(i-1)+1:4*i,4*(j-1)+1:4*j);
        
        redblock_H=(1/4)*H_N*redblock;
        greenblock_H=(1/4)*H_N*greenblock;
        blueblock_H=(1/4)*H_N*blueblock;
        
        avg_12_red=(redblock_H(1,1)+redblock_H(1,2))/2;
        avg_34_red=(redblock_H(1,3)+redblock_H(1,4))/2;
        
        avg_12_green=(greenblock_H(1,1)+greenblock_H(1,2))/2;
        avg_34_green=(greenblock_H(1,3)+greenblock_H(1,4))/2;
        
        avg_12_blue=(blueblock_H(1,1)+blueblock_H(1,2))/2;
        avg_34_blue=(blueblock_H(1,3)+blueblock_H(1,4))/2;
        
        
        T=5;d=5;
        redblock_H_E=redblock_H;
        greenblock_H_E=greenblock_H;
        blueblock_H_E=blueblock_H;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RED EMBEDDING
         
        if(watermark_scramble_red((i-1)*256+(2*j)-1)==1)
            if(redblock_H(1,1)-redblock_H(1,2)<=d)
                redblock_H_E(1,1)=avg_12_red+(T/2);
            end
        end
        if(watermark_scramble_red((i-1)*256+(2*j)-1)==0)
            if(redblock_H(1,2)-redblock_H(1,1)<=d)
                redblock_H_E(1,1)=avg_12_red-(T/2);
            end
        end
        
           
        if(watermark_scramble_red((i-1)*256+(2*j)-1)==0)
            if(redblock_H(1,2)-redblock_H(1,1)<=d)
                redblock_H_E(1,2)=avg_12_red+(T/2);
            end
        end
        if(watermark_scramble_red((i-1)*256+(2*j)-1)==1)
            if(redblock_H(1,1)-redblock_H(1,2)<=d)
                redblock_H_E(1,2)=avg_12_red-(T/2);
            end
        end
        
        
          
        if(watermark_scramble_red((i-1)*256+(2*j))==1)
            if(redblock_H(1,3)-redblock_H(1,4)<=d)
                redblock_H_E(1,3)=avg_34_red+(T/2);
            end
        end
        if(watermark_scramble_red((i-1)*256+(2*j))==0)
            if(redblock_H(1,4)-redblock_H(1,3)<=d)
                redblock_H_E(1,3)=avg_34_red-(T/2);
            end
        end
        
           
        if(watermark_scramble_red((i-1)*256+(2*j))==0)
            if(redblock_H(1,4)-redblock_H(1,3)<=d)
                redblock_H_E(1,4)=avg_34_red+(T/2);
            end
        end
        if(watermark_scramble_red((i-1)*256+(2*j))==1)
            if(redblock_H(1,3)-redblock_H(1,4)<=d)
                redblock_H_E(1,4)=avg_34_red-(T/2);
            end
        end
        
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GREEN EMBEDDING
         
        if(watermark_scramble_green((i-1)*256+(2*j)-1)==1)
            if(greenblock_H(1,1)-greenblock_H(1,2)<=d)
                greenblock_H_E(1,1)=avg_12_green+(T/2);
            end
        end
        if(watermark_scramble_green((i-1)*256+(2*j)-1)==0)
            if(greenblock_H(1,2)-greenblock_H(1,1)<=d)
                greenblock_H_E(1,1)=avg_12_green-(T/2);
            end
        end
        
           
        if(watermark_scramble_green((i-1)*256+(2*j)-1)==0)
            if(greenblock_H(1,2)-greenblock_H(1,1)<=d)
                greenblock_H_E(1,2)=avg_12_green+(T/2);
            end
        end
        if(watermark_scramble_green((i-1)*256+(2*j)-1)==1)
            if(greenblock_H(1,1)-greenblock_H(1,2)<=d)
                greenblock_H_E(1,2)=avg_12_green-(T/2);
            end
        end
        
        
          
        if(watermark_scramble_green((i-1)*256+(2*j))==1)
            if(greenblock_H(1,3)-greenblock_H(1,4)<=d)
                greenblock_H_E(1,3)=avg_34_green+(T/2);
            end
        end
        if(watermark_scramble_green((i-1)*256+(2*j))==0)
            if(greenblock_H(1,4)-greenblock_H(1,3)<=d)
                greenblock_H_E(1,3)=avg_34_green-(T/2);
            end
        end
        
           
        if(watermark_scramble_green((i-1)*256+(2*j))==0)
            if(greenblock_H(1,4)-greenblock_H(1,3)<=d)
                greenblock_H_E(1,4)=avg_34_green+(T/2);
            end
        end
        if(watermark_scramble_green((i-1)*256+(2*j))==1)
            if(greenblock_H(1,3)-greenblock_H(1,4)<=d)
                greenblock_H_E(1,4)=avg_34_green-(T/2);
            end
        end
        
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           
           
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLUE EMBEDDING
         
        if(watermark_scramble_blue((i-1)*256+(2*j)-1)==1)
            if(blueblock_H(1,1)-blueblock_H(1,2)<=d)
                blueblock_H_E(1,1)=avg_12_blue+(T/2);
            end
        end
        if(watermark_scramble_blue((i-1)*256+(2*j)-1)==0)
            if(blueblock_H(1,2)-blueblock_H(1,1)<=d)
                blueblock_H_E(1,1)=avg_12_blue-(T/2);
            end
        end
        
           
        if(watermark_scramble_blue((i-1)*256+(2*j)-1)==0)
            if(blueblock_H(1,2)-blueblock_H(1,1)<=d)
                blueblock_H_E(1,2)=avg_12_blue+(T/2);
            end
        end
        if(watermark_scramble_blue((i-1)*256+(2*j)-1)==1)
            if(blueblock_H(1,1)-blueblock_H(1,2)<=d)
                blueblock_H_E(1,2)=avg_12_blue-(T/2);
            end
        end
        
        
          
        if(watermark_scramble_blue((i-1)*256+(2*j))==1)
            if(blueblock_H(1,3)-blueblock_H(1,4)<=d)
                blueblock_H_E(1,3)=avg_34_blue+(T/2);
            end
        end
        if(watermark_scramble_blue((i-1)*256+(2*j))==0)
            if(blueblock_H(1,4)-blueblock_H(1,3)<=d)
                blueblock_H_E(1,3)=avg_34_blue-(T/2);
            end
        end
        
           
        if(watermark_scramble_blue((i-1)*256+(2*j))==0)
            if(blueblock_H(1,4)-blueblock_H(1,3)<=d)
                blueblock_H_E(1,4)=avg_34_blue+(T/2);
            end
        end
        if(watermark_scramble_blue((i-1)*256+(2*j))==1)
            if(blueblock_H(1,3)-blueblock_H(1,4)<=d)
                blueblock_H_E(1,4)=avg_34_blue-(T/2);
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
        
        
        redblock_E=H_N*redblock_H_E;
        greenblock_E=H_N*greenblock_H_E;
        blueblock_E=H_N*blueblock_H_E;
        
        watermarked(4*(i-1)+1:4*i,4*(j-1)+1:4*j,1)=redblock_E;
        watermarked(4*(i-1)+1:4*i,4*(j-1)+1:4*j,2)=greenblock_E;
        watermarked(4*(i-1)+1:4*i,4*(j-1)+1:4*j,3)=blueblock_E;
    end
end
  watermarked=uint8(watermarked);
  
%   watermarked=imnoise(watermarked,'salt & pepper',0.01);
%   
%   watermarked=imnoise(watermarked,'gaussian',0,0.001);
%   
%   watermarked=imnoise(watermarked,'speckle',0.001);
%   

%   watermarked(:,:,1)=medfilt2(watermarked(:,:,1),[2,1]);
%   watermarked(:,:,2)=medfilt2(watermarked(:,:,2),[2,1]);
%   watermarked(:,:,3)=medfilt2(watermarked(:,:,3),[2,1]);
%   

%   watermarked=imrotate(watermarked,5);watermarked=imresize(watermarked,[512,512]);
%   
%   watermarked=imresize(watermarked,[1024,1024]);cover=imresize(cover,[1024,1024]);
%   
%   imwrite(uint8(watermarked),'watermarkd.jpg');watermarked=imread('watermarkd.jpg');
 
  watermarked=double(watermarked);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRACTION %%%%%%%%%%%%%%%%%%%%%%%

for i = 1:128
    for j =1:128
        redblock_E=watermarked(4*(i-1)+1:4*i,4*(j-1)+1:4*j,1);
        greenblock_E=watermarked(4*(i-1)+1:4*i,4*(j-1)+1:4*j,2);
        blueblock_E=watermarked(4*(i-1)+1:4*i,4*(j-1)+1:4*j,3);
        redblock_E_EX=(1/4)*H_N*redblock_E;
        if(redblock_E_EX(1,1)<=redblock_E_EX(1,2))
            watermark_scramble_red_extract((i-1)*256+(2*j)-1)=0;
        else
            watermark_scramble_red_extract((i-1)*256+(2*j)-1)=1;                     % RED EXTRACTION
        end
        
        if(redblock_E_EX(1,3)<=redblock_E_EX(1,4))
            watermark_scramble_red_extract((i-1)*256+(2*j))=0;
        else
            watermark_scramble_red_extract((i-1)*256+(2*j))=1;
        end
        
        
        greenblock_E_EX=(1/4)*H_N*greenblock_E;
        if(greenblock_E_EX(1,1)<=greenblock_E_EX(1,2))
            watermark_scramble_green_extract((i-1)*256+(2*j)-1)=0;
        else
            watermark_scramble_green_extract((i-1)*256+(2*j)-1)=1;                     % GREEN EXTRACTION
        end
        
        if(greenblock_E_EX(1,3)<=greenblock_E_EX(1,4))
            watermark_scramble_green_extract((i-1)*256+(2*j))=0;
        else
            watermark_scramble_green_extract((i-1)*256+(2*j))=1;
        end
          
          blueblock_E_EX=(1/4)*H_N*blueblock_E;
        if(blueblock_E_EX(1,1)<=blueblock_E_EX(1,2))
            watermark_scramble_blue_extract((i-1)*256+(2*j)-1)=0;
        else
            watermark_scramble_blue_extract((i-1)*256+(2*j)-1)=1;                     % BLUE EXTRACTION
        end
        
        if(blueblock_E_EX(1,3)<=blueblock_E_EX(1,4))
            watermark_scramble_blue_extract((i-1)*256+(2*j))=0;
        else
            watermark_scramble_blue_extract((i-1)*256+(2*j))=1;
        end
        
    end
end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%

watermark_scramble_red_extract=reshape(watermark_scramble_red_extract,[64*64,8]);
watermark_scramble_green_extract=reshape(watermark_scramble_green_extract,[64*64,8]);
watermark_scramble_blue_extract=reshape(watermark_scramble_blue_extract,[64*64,8]);

watermark_scramble_red_extract=bi2de(watermark_scramble_red_extract);
watermark_scramble_green_extract=bi2de(watermark_scramble_green_extract);
watermark_scramble_blue_extract=bi2de(watermark_scramble_blue_extract);

watermark_scramble_red_extract=reshape(watermark_scramble_red_extract,[64,64]);
watermark_scramble_green_extract=reshape(watermark_scramble_green_extract,[64,64]);
watermark_scramble_blue_extract=reshape(watermark_scramble_blue_extract,[64,64]);

for i=0:63
    for j=0:63
        coordinator=([2,-1;-1,1]*[i;j])+[64;64];
        m=rem(coordinator(1),64);
        n=rem(coordinator(2),64);
        watermark_rescramble_red(m+1,n+1)=watermark_scramble_red_extract(i+1,j+1);
    end
end


for i=0:63
    for j=0:63
        coordinator=([2,-1;-1,1]*[i;j])+[64;64];
        m=rem(coordinator(1),64);
        n=rem(coordinator(2),64);
        watermark_rescramble_blue(m+1,n+1)=watermark_scramble_blue_extract(i+1,j+1);
    end
end


for i=0:63
    for j=0:63
        coordinator=([2,-1;-1,1]*[i;j])+[64;64];
        m=rem(coordinator(1),64);
        n=rem(coordinator(2),64);
        watermark_rescramble_green(m+1,n+1)=watermark_scramble_green_extract(i+1,j+1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




watermark_rescramble(:,:,1)=watermark_rescramble_red;
watermark_rescramble(:,:,2)=watermark_rescramble_green;
watermark_rescramble(:,:,3)=watermark_rescramble_blue;

watermark_scramble_red=reshape(watermark_scramble_red,[64*64,8]);
watermark_scramble_green=reshape(watermark_scramble_green,[64*64,8]);
watermark_scramble_blue=reshape(watermark_scramble_blue,[64*64,8]);

watermark_scramble_red=bi2de(watermark_scramble_red);
watermark_scramble_green=bi2de(watermark_scramble_green);
watermark_scramble_blue=bi2de(watermark_scramble_blue);

watermark_scramble_red=reshape(watermark_scramble_red,[64,64]);
watermark_scramble_green=reshape(watermark_scramble_green,[64,64]);
watermark_scramble_blue=reshape(watermark_scramble_blue,[64,64]);

watermark_scramble(:,:,1)=watermark_scramble_red;
watermark_scramble(:,:,2)=watermark_scramble_green;
watermark_scramble(:,:,3)=watermark_scramble_blue;


MSE=mean(mean(mean((cover-watermarked).^2)))
PSNR=10*log10((255^2)/MSE)
EMBEDDING_CAPACITY=(size(cover,1)*size(cover,2))/(size(watermark,1)*size(watermark,2))
NC=0;
for i=1:64
    for j=1:64
        if(watermark(i,j)==watermark_rescramble(i,j))
        NC=NC+1;
        end
    end
end
normalizedcoefficient=NC/(64*64)
structuralsimilarity=ssim(cover,watermarked)


imshow(uint8(cover));
figure, imshow(uint8(watermarked));
figure,imshow(uint8(watermark));
figure,imshow(uint8(watermark_scramble));
figure,imshow(uint8(watermark_rescramble));