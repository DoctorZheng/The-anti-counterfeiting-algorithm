function  result = RALBGC(I)
%I=[44 44 44;44 12 5; 5 5 5];
[m,n,h] = size(I);
if h==3
    I =  rgb2gray(I);
end
I=double(I);
lbpIu = double(zeros([m-2 n-2]));
lbpId = double(zeros([m-2 n-2]));
t=0;
medB=mean(mean(I));
for i = 2:m-1
    for j = 2:n-1
        medS=mean(mean(I(i-1:i+1,j-1:j+1)));
        neighbor01 = ([medS I(i,j-1) I(i+1,j-1) I(i+1,j) I(i+1,j+1) I(i,j+1) I(i-1,j+1) I(i-1,j) I(i-1,j-1)] - [I(i,j) I(i+1,j-1) I(i+1,j) I(i+1,j+1) I(i,j+1) I(i-1,j+1) I(i-1,j) I(i-1,j-1) I(i,j-1)])>t;
        neighbor02 = ([medB I(i+1,j) I(i+1,j+1) I(i,j+1) I(i-1,j+1) I(i-1,j) I(i-1,j-1) I(i,j-1) I(i+1,j-1)] - [I(i,j) I(i+1,j-1) I(i+1,j) I(i+1,j+1) I(i,j+1) I(i-1,j+1) I(i-1,j) I(i-1,j-1) I(i,j-1)])>t;
        neighbor_a = (neighbor01 + neighbor02) > 1;
        neighbor03 = ([medS I(i,j-1) I(i+1,j-1) I(i+1,j) I(i+1,j+1) I(i,j+1) I(i-1,j+1) I(i-1,j) I(i-1,j-1)] - [I(i,j) I(i+1,j-1) I(i+1,j) I(i+1,j+1) I(i,j+1) I(i-1,j+1) I(i-1,j) I(i-1,j-1) I(i,j-1)])<=t;
        neighbor04 = ([medB I(i+1,j) I(i+1,j+1) I(i,j+1) I(i-1,j+1) I(i-1,j) I(i-1,j-1) I(i,j-1) I(i+1,j-1)] - [I(i,j) I(i+1,j-1) I(i+1,j) I(i+1,j+1) I(i,j+1) I(i-1,j+1) I(i-1,j) I(i-1,j-1) I(i,j-1)])<=t;
        neighbor_r = (neighbor03 + neighbor04) > 1;
        pixel_u = 0;pixel_d = 0;
        for k = 1:9
            pixel_u = pixel_u + neighbor_a(1,k) * bitshift(1,9-k);
            pixel_d = pixel_d + neighbor_r(1,k) * bitshift(1,9-k);
        end
        lbpIu(i-1,j-1) = double(pixel_u);
        lbpId(i-1,j-1) = double(pixel_d);
    end
end

% Return with LBP histogram
bins= 256;%8*(8-1) + 3;  
% result(1,:)=hist(lbpIu(:),0:(bins-1));%/ numel(lbpIu);
% result(2,:)=hist(lbpId(:),0:(bins-1));%/ numel(lbpIu);
% result(1,:)=result(1,:)/sum(result(1,:));
% result(2,:)=result(2,:)/sum(result(2,:));
result1=hist(lbpIu(:),0:(bins-1));%/ numel(lbpIu);
result2=hist(lbpId(:),0:(bins-1));%/ numel(lbpIu);
result=[result1,result2];
result=result/sum(result);