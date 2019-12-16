clear;
clc;
tic;
%I=imread('C:\Users\96446\Desktop\rgb\1305031910.865025.png');
%I=imread('C:\Users\96446\Desktop\rgb\1305031912.197357.png');
I=imread('D:\论文\边线提取\pic\6.png');
%I=imread(rgb);
figure;
imshow(I);
I1=rgb2gray(I);
%h=fspecial('gaussian',5);
I2=edge(I1,'canny');
figure;
label=255*I2;
image(label);
xlim([1 size(label,2)]);
ylim([1 size(label,1)]);
toc;
%[H,T,R] = hough(I2);
%PEAKS = houghpeaks(H,100);
%lines=houghlines(I2,T,R,PEAKS);
%for k=1:length(lines)
    %xy=[lines(k).point1;lines(k).point2];
    %plot(xy(:,1),xy(:,2),'Linewidth',3);
%end
