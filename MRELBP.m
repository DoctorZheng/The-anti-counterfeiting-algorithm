function  result = MRELBP(J)
[m,n,h] = size(J);
if h==3
    J = rgb2gray(J);
end
I=medfilt2(J,[5,5]); 
Ismall=medfilt2(J,[3,3]); 
ICrop = I(2:m-1,2:n-1);
ImeanCrop = mean(mean(ICrop));
ImageCenterCI = (I - ImeanCrop) > 0;
bins= 2;%8*(8-1) + 3;  
resultCI=hist(ImageCenterCI(:),0:(bins-1));

ImageCenterNI = double(zeros([m-2 n-2]));
for i = 2:m-1
    for j = 2:n-1
        ICropNI = I(i-1:i+1,j-1:j+1);        
        ImeanCropNI=mean(mean(ICropNI));
        CenterNI = ICropNI - ImeanCropNI;
        pixel = 0;
        NI=[CenterNI(2,3) CenterNI(1,3) CenterNI(1,2) CenterNI(1,1) CenterNI(2,2-1) CenterNI(2+1,2-1) CenterNI(2+1,2) CenterNI(2+1,2+1)] > 0;%ÄæÊ±Õë
        for k = 1:8
            pixel = pixel + NI(1,k) * bitshift(1,k-1);
        end
        ImageCenterNI(i-1,j-1) = double(pixel);
        
        ICropNIsmall = Ismall(i-1:i+1,j-1:j+1);
        CenterRD = ICropNI - ICropNIsmall;
        pixel = 0;
        RD=[CenterRD(2,2+1) CenterRD(2-1,2+1) CenterRD(2-1,2) CenterRD(2-1,2-1) CenterRD(2,2-1) CenterRD(2+1,2-1) CenterRD(2+1,2) CenterRD(2+1,2+1)] > 0;%ÄæÊ±Õë
        for k = 1:8
            pixel = pixel + RD(1,k) * bitshift(1,k-1);
        end
        ImageCenterRD(i-1,j-1) = double(pixel);        
    end
end

% Return with LBP histogram
bins= 256;%8*(8-1) + 3;  
resultNI=hist(ImageCenterNI(:),0:(bins-1));
resultRD=hist(ImageCenterRD(:),0:(bins-1));

result=[resultCI,resultNI,resultRD];
result=result/sum(result);

