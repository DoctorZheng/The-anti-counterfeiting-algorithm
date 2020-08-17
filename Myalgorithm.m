close all;
clear all;

file_path =  '.\1\';% 图像文件夹路径，
cc=colormap(lines(100));
img_path_list = dir(strcat(file_path,'*.bmp'));%获取该文件夹中所有jpg格式的图像。
img_num = length(img_path_list);%获取图像总数量，
NUMD=zeros(img_num,10);
for ii = 1:10 %逐一读取图像，    
     I=imread(strcat(num2str(ii),'.bmp'));
     file_path =  strcat(strcat('.\',num2str(ii)),'\');
    %SVDis=zeros(size(STATS0,1),img_num);
    if img_num > 0 %有满足条件的图像，
        for jj = 1:img_num %逐一读取图像，
            image_name = img_path_list(jj).name;% 图像名，
            I1 =  imread(strcat(file_path,image_name));
            
            thOTSU = graythresh(I);
            imagBWO = im2bw(I, thOTSU);
            thOTSU1 = graythresh(I1);
            imagBWO1 = im2bw(I1, thOTSU1);
            [H W]=size(I);

            %差异区域提取
            imagBWO2 = imagBWO1;
            I0=double(I);
            I2=double(I1);
            choice=1;
            switch choice
                case 1
                    for i = 1 : H
                        for j = 1 : W
                            if (imagBWO(i,j)==1&&imagBWO1(i,j)==0)||(imagBWO(i,j)==0&&imagBWO1(i,j)==1)
                                tcc = abs(I2(i,j)-thOTSU1*255);
                                if tcc<30  %I1(i,j)<thOTSU1*255+k * sigma
                                    imagBWO2(i,j)=~imagBWO2(i,j);
                                    %imagBWO2(i,j)=255;
                                end
                            end
                        end
                    end
                case 2
                    b=3;
                    k = 0.34;%-0.2;
                    % Mean value
                    window=[15 15];
                    padding = 'replicate';
                    mean1 = averagefilter(I2, window, padding);
                    % Standard deviation
                    meanSquare = averagefilter(I2.^2, window, padding);
                    deviation = (meanSquare - mean1.^2).^0.5;
                    % Sauvola
                    R = max(deviation(:));
                    threshold = mean1.*(power((deviation./256+2)*0.4,2));%mean1.*(1 + k * (deviation / R-1));
                    %threshold = mean1.*(1 + 0.05 * ((I2-mean1)/(1-I2+mean1)-1));
                    for i = 1 : H
                        for j = 1 : W
                            if (imagBWO(i,j)==1&&imagBWO1(i,j)==0)||(imagBWO(i,j)==0&&imagBWO1(i,j)==1)
                                % upR =  i-floor(b/2-1/2);
                                % dnR =  i+floor(b/2);
                                % lfC =  j-floor(b/2-1/2);
                                % rtC =  j+floor(b/2);
                                % m_ij = mean(mean(I1(upR : dnR, lfC : rtC)));
                                % sigma_squared = double(I1(upR : dnR, lfC : rtC)) - repmat(m_ij, (dnR-upR+1), (rtC-lfC+1));
                                % sigma_squared = sigma_squared .^ 2;
                                % sigma_squared = mean(mean(sigma_squared));
                                % sigma = sqrt(sigma_squared);
                                % t_ij = m_ij.*(1 + k * (sigma / R-1));
                                if I2(i,j)<threshold(i,j)
                                    imagBWO2(i,j)=0;
                                else
                                    imagBWO2(i,j)=255;
                                end
                            end
                        end
                    end
                case 3
                    k = 0.3;%-0.2;
                    % Mean value
                    window=[15 15];
                    padding = 'replicate';
                    mean0 = averagefilter(I0, window, padding);
                    % Standard deviation
                    meanSquare0 = averagefilter(I0.^2, window, padding);
                    deviation0 = (meanSquare0 - mean0.^2).^0.5;
                    % Sauvola
                    R0 = max(deviation0(:));
                    threshold0 = mean0.*(1 + k * (deviation0 / R0-1));        
                    mean1 = averagefilter(I2, window, padding);
                    % Standard deviation
                    meanSquare1 = averagefilter(I2.^2, window, padding);
                    deviation1 = (meanSquare1 - mean1.^2).^0.5;
                    % Sauvola
                    R1 = max(deviation1(:));
                    threshold1 = mean1.*(1 + k * (deviation1 / R1-1));

                    for i = 1 : H
                        for j = 1 : W
                            if imagBWO(i,j)~= imagBWO1(i,j)%(imagBWO(i,j)==1&&imagBWO1(i,j)==0)||(imagBWO(i,j)==0&&imagBWO1(i,j)==1)
                                D0=I0(i,j)-mean0(i,j)*(power((deviation0(i,j)/256+1)*0.5,2));%*k-(1-k)*(deviation0(i,j));
                                D1=I2(i,j)-mean1(i,j)*(power((deviation0(i,j)/256+1)*0.5,2));%*k-(1-k)*(deviation1(i,j));
                                if D0/D1>0%(I0(i,j)-threshold0(i,j))/ (I2(i,j)-threshold1(i,j))>0  %D0/D1>0&&abs(D0-D1)<100
                                    imagBWO2(i,j)=~imagBWO2(i,j);
                                else
                                    bb=0;
                                end
                            end
                        end
                    end      
                case 4
                    b=3;
                    for i = 1+b : H-b
                        for j = 1+b: W-b
                            if (imagBWO(i,j)==1&&imagBWO1(i,j)==0)||(imagBWO(i,j)==0&&imagBWO1(i,j)==1)
                                upR =  i-floor(b);
                                dnR =  i+floor(b);
                                lfC =  j-floor(b);
                                rtC =  j+floor(b);
                                C1=I(upR : dnR, lfC : rtC);
                                C2=I1(upR : dnR, lfC : rtC);
                                thOTSU1 = graythresh(C1);
                                thOTSU2 = graythresh(C2);
                                if (I(i,j)<thOTSU1&&I1(i,j)<thOTSU2)||(I(i,j)>thOTSU1&&I1(i,j)>thOTSU2)
                                    imagBWO2(i,j)=imagBWO(i,j);
                                end
                            end
                        end
                    end
                otherwise
                    display('Wrong Choice!');
            end

            %删除小面积对象
            imagBWO = bwareaopen(imagBWO,6);
            imagBWO2 = bwareaopen(imagBWO2,6);
t1=cputime;
            %骨架提取
            bw0=bwmorph(imagBWO,'skel',Inf);
            bw2=bwmorph(bw0,'spur',1);

            %骨架梯度   
            [Fx,Fy]=gradient(double(bw2));
            %g1=sqrt(Fx.^2+Fy.^2);
            Iphase=atan(Fy./Fx).*180./pi;
            Iphase(Iphase(:,:)==-90) = 90;

            %骨架法线方向
            [bonePointRows, bonePointCols] = find(bw2);
            bw3=double(bw2)*nan;
            bw4=double(bw2)*nan;
            for i=1:size(bonePointRows)
                %使用3*3的窗口对直线两边的梯度进行统计分析，获取最大次数的梯度作为该点法线   
                A=Iphase( max(1,bonePointRows(i)-1):min(H,bonePointRows(i)+1),max(1,bonePointCols(i)-1):min(W,bonePointCols(i)+1));         
                B=reshape(A.',1,numel(A));
                %出现次数最多
                %含有nan值矩阵均值
                notNanValues=B(~isnan(B));
                X1 = unique(notNanValues);
                if length(X1)==1
                    [M1 N1]=hist(notNanValues,1);
                else
                    [M1 N1]=hist(notNanValues,X1);
                end
                maxcount = max(M1);   
                T2=find(M1(:)==maxcount);
                result=X1(T2);
                meanvalue=sum(notNanValues)./length(notNanValues);
                [mi,q]=min(abs(result-meanvalue));  
                angle=result(q);
                bw3(bonePointRows(i),bonePointCols(i))=angle;
                %计算该点骨骼宽度
                maxStrokeWidth=50;
                step0=1;
                while step0 < maxStrokeWidth
                    nextR = round(bonePointRows(i) + sin(angle/180*pi) * step0);
                    nextC = round(bonePointCols(i) + cos(angle/180*pi) * step0);
                    step0 = step0 + 1;
                    if nextR < 1 | nextC < 1 | nextR > H | nextC > W
                        break
                    end  
                    if imagBWO(nextR,nextC)==0
                        break
                    end    
                end
                step1=1;
                while step1 < maxStrokeWidth
                    nextR = round(bonePointRows(i) + sin(angle/180*pi) * step1);
                    nextC = round(bonePointCols(i) + cos(angle/180*pi) * step1);
                    step1 = step1 + 1;
                    if nextR < 1 | nextC < 1 | nextR > H | nextC > W
                        break
                    end  
                    if imagBWO(nextR,nextC)==0
                        break
                    end    
                end
                bw4(bonePointRows(i),bonePointCols(i))=step0+step1-1;
            end
e1=cputime-t1;
            %剔除掉单个3以下微小宽度点。
            bw5=imagBWO*0;
            for i=1:size(bonePointRows)
                if bw4(bonePointRows(i),bonePointCols(i))<4        
                    A=bw4(max(1,bonePointRows(i)-1):min(H,bonePointRows(i)+1),max(1,bonePointCols(i)-1):min(W,bonePointCols(i)+1));         
                    A(2,2)=nan;
                    [m,n]=find(~isnan(A));
                    if numel(A)~=9|numel(m)<2|(numel(m)==2&((m(1)==m(2)&abs(n(1)-n(2))==1)|(n(1)==n(2)&abs(m(1)-m(2))==1)))
                       bw4(bonePointRows(i),bonePointCols(i))= 4;%mean
                    else  
                        bw5(bonePointRows(i),bonePointCols(i))= 1;
                    end
                end
            end

            %细微区域连通域，最好用24连通即5*5
            %bw5=logical(bw5);
             L0 = bwlabel(bw5,8);
             STATS0 = regionprops(L0,'all');%STATS中含有所有连通域的properations

            %样本连通域
            L = bwlabel(imagBWO);

            %待测连通域
            L1 = bwlabel(imagBWO2);

            %判断差异
            NUMDIFF=0;
            for i = 1 : size(STATS0, 1)
                  pointlist = STATS0(i).PixelList;
                  boundary0 = STATS0(i).BoundingBox;
                  boundary0(1,1)=max(0,boundary0(1,1)-2);
                  boundary0(1,2)=max(0,boundary0(1,2)-2);
                  boundary0(1,3)=min(W-boundary0(1,1),boundary0(1,3)+4);
                  boundary0(1,4)=min(H-boundary0(1,2),boundary0(1,4)+4);
                  Image1=imcrop(L,boundary0);
                  Image2=imcrop(L1,boundary0);
                  statisl = tabulate(reshape(Image1.',1,numel(Image1)));
                  statis2 = tabulate(reshape(Image2.',1,numel(Image2)));
                  num1=sum(statisl(:,1)~=0);
                  num2=sum(statis2(:,1)~=0);
                  if(num1~=num2)
                      NUMDIFF=NUMDIFF+1;
                  end
            end
            NUMD(jj,ii)=NUMDIFF;
        end
    end
end


