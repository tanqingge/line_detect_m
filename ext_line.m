%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
num_segments=1;
MAX_Segments=10000;
ni=size(label,1);
nj=size(label,2);
widthMin=1;
widthMax=5;
connect_th=3;
max_gap=10;
min_length=10;
state=0;
width=0;
width_1=0;
%{segment_init:}%
for i=1:MAX_Segments
segments(i).state=0;
segments(i).gap=0;
segments(i).count=0;
segments(i).xy=0;
segments(i).xx=0;
segments(i).x=0;
segments(i).y=0;
segments(i).updated=0;  
segments(i).grad=0;
segments(i).max_width=0;
end
%Scan for vertical field line pixels:
for j=1:nj
    for i=1:ni
        label_single=label(i,j);
           switch state
               case 0
               if label_single==0
                   state =1;
                   width_1=0;
               else
                   width_1=0;
               end
               case 1
               if label_single==255
                   state=2;
                   width=1;
                   width_1=0;
               elseif label~=0
                   state=0;
                   width_1=0;
               else 
                   width_1=0;
               end
               case 2
                   if label_single==0
                       state=1;
                       width_1=width;
                   elseif label~=255
                       state=0;
                       width_1=0;
                   else
                       width=width+1;
                       width_1=0; 
                   end
           end
        if ((width_1>=widthMin)&&(width_1<=widthMax))
            iline=i-(width_1+1)/2;
            seg_updated=0;
            for k=1:num_segments
                if segments(k).state==1
                    xProj = double(segments(k).xMean + segments(k).invgrad*(j-segments(k).yMean));
                    xErr = double(iline-xProj);
                    if xErr<0
                        xErr=-xErr;
                    end
                    if xErr<connect_th
                        %updateStat(k,i,j,max_gap)%
                        segments(k).xMean=(segments(k).count*segments(k).xMean+iline)/(segments(k).count+1);
                        segments(k).yMean=(segments(k).count*segments(k).yMean+j)/(segments(k).count+1);
                        segments(k).x=segments(k).x+iline;
                        segments(k).y=segments(k).y+j;
                        segments(k).xy=segments(k).xy+iline*j;
                        segments(k).xx=segments(k).xx+iline*iline;
                        segments(k).yy=segments(k).yy+j*j;
                        segments(k).gap=segments(k).gap+1;
                        if segments(k).gap>(max_gap+1)
                            segments(k).gap=max_gap+1;
                        end
                        segments(k).count=segments(k).count+1;
                        segments(k).invgrad=double((segments(k).xy- segments(k).x*segments(k).y/segments(k).count)/(segments(k).yy-segments(k).y*segments(k).y/segments(k).count));
                        if (segments(k).invgrad>1.0)||(segments(k).invgrad<-1.0)
                            segments(k).state=2; 
                            segments(k).count=0;
                        end
                        segments(k).updated=1;
                        segments(k).x1=iline;
                        segments(k).y1=j;
                        seg_updated=seg_updated+1;
                        if (width > segments(k).max_width)
                            segments(k).max_width = width;
                        end
                    end
                end
            end
            if (seg_updated==0)&&(num_segments<MAX_Segments)
                %initStat%
                segments(num_segments).state=1;
                segments(num_segments).count=1;
                segments(num_segments).gap=1;
                segments(num_segments).grad=0;
                segments(num_segments).x0=iline;
                segments(num_segments).y0=j;
                segments(num_segments).x=iline;
                segments(num_segments).y=j;
                segments(num_segments).xx=iline*iline;
                segments(num_segments).xy=iline*j;
                segments(num_segments).yy=j*j;
                segments(num_segments).xMean=iline;
                segments(num_segments).yMean=j;
                segments(num_segments).updated=1;
                segments(num_segments).invgrad=0;
                segments(num_segments).length=0;
                segments(num_segments).max_width=10;   
                num_segments=num_segments+1;
            end
        end
    end
 %Segments refresh%
    for i_1=1:num_segments
        if (segments(i_1).state==1)&&(segments(i_1).updated==0)
            if segments(i_1).gap>0
                segments(i_1).gap=segments(i_1).gap-1;
            else
                segments(i_1).state=2;
            end
        end
        segments(i_1).updated=0;
    end
end
 %Segments terminate%
for i=1:num_segments
    segments(i).state=2;
end
%Scan for horizontal field line pixels:
for i=1:ni
    for j=1:nj
    label_single=label(i,j);
%Linestate%
           switch state
               case 0
               if label_single==0
                   state =1;
                   width_1=0;
               else
                   width_1=0;
               end
               case 1
               if label_single==255
                   state=2;
                   width=1;
                   width_1=0;
               elseif label~=0
                   state=0;
                   width_1=0;
               else
                   width_1=0;
               end
               case 2'
                   if label_single==0
                       state=1;
                       width_1=width;
                   elseif label~=255
                       state=0;
                       width_1=0;
                   else
                       width=width+1;
                       width_1=0;
                   end
           end
    if ((width_1>=widthMin)&&(width_1<=widthMax))
        jline=j-(width_1+1)/2;
        seg_updated=0;
        for k=1:num_segments
            if segments(k).state==1
                yProj=double(segments(k).yMean+segments(k).grad*(i-segments(k).xMean));
                yErr=double(jline-yProj);
                if yErr<0
                    yErr=-yErr;
                end
                if yErr<connect_th
                    %updateStat(k,i,j,max_gap):
                    segments(k).xMean=(segments(k).count*segments(k).xMean+i)/(segments(k).count+1);
                    segments(k).yMean=(segments(k).count*segments(k).yMean+jline)/(segments(k).count+1);
                    segments(k).x=segments(k).x+i;
                    segments(k).y=segments(k).y+jline;
                    segments(k).xy=segments(k).xy+i*jline;
                    segments(k).xx=segments(k).xx+i*i;
                    segments(k).yy=segments(k).yy+jline*jline;
                    segments(k).gap=segments(k).gap+1;
                    if segments(k).gap>(max_gap+1)
                        segments(k).gap=max_gap+1;
                    end
                    segments(k).count=segments(k).count+1;
                    segments(k).grad=double((segments(k).xy-segments(k).x*segments(k).y/segments(k).count)/(segments(k).xx-segments(k).x*segments(k).x/segments(k).count));
                    if (segments(k).grad>1)||(segments(k).grad<-1)
                        segments(k).state=2;
                        segments(k).count=0;
                    end
                    segments(k).updated=1;
                    segments(k).x1=i;
                    segments(k).y1=jline;
                    if width>segments(k).max_width
                        segments(k).max_width=width;
                    end
                    seg_updated=seg_updated+1;
                end
            end
        end
        if (seg_updated==0)&&(num_segments<MAX_Segments)
            segments(num_segments).state=1;
            segments(num_segments).count=1;
            segments(num_segments).gap=1;
            segments(num_segments).grad=0;
            segments(num_segments).x0=i;
            segments(num_segments).y0=jline;
            segments(num_segments).x=i;
            segments(num_segments).y=jline;
            segments(num_segments).xx=i*i;
            segments(num_segments).xy=i*jline;
            segments(num_segments).yy=jline*jline;
            segments(num_segments).xMean=i;
            segments(num_segments).yMean=jline;
            segments(num_segments).updated=1;
            segments(num_segments).invgrad=0;
            segments(num_segments).length=0;
            segments(num_segments).max_width=10;
            num_segments=num_segments+1;
        end
    end
    end
%segment_refresh:
    for i_2=1:num_segments
    if (segments(i_2).state==1)&&(segments(i_2).updated==0)
        if segments(i_2).gap>0
            segments(i_2).gap=segments(i_2).gap-1;
        else
            segments(i_2).state=2;
        end
    end
    segments(i_2).updated=0;
    end
end
valid_segments=0;
figure;
for r=1:num_segments
    segments(r).dx=segments(r).x1-segments(r).x0;
    segments(r).dy=segments(r).y1-segments(r).y0;
    segments(r).length = sqrt(segments(r).dx*segments(r).dx+segments(r).dy*segments(r).dy);
    if (segments(r).count>min_length)
        line([segments(r).y0 segments(r).y1],[segments(r).x0 segments(r).x1]);
        valid_segments=valid_segments+1;
        available_segments(valid_segments).x0=segments(r).x0;
        available_segments(valid_segments).y0=segments(r).y0;
        available_segments(valid_segments).x1=segments(r).x1;
        available_segments(valid_segments).y1=segments(r).y1;
        available_segments(valid_segments).grad=segments(r).grad;
        available_segments(valid_segments).invgrad=segments(r).invgrad;
        
        hold on;
    end
end
xlim([1 size(label,2)]);
ylim([1 size(label,1)]);
length=5;
X=0:length:nj;
Y=0:length:ni;
%ndex_X=232;
%index_Y=148;    
%M=meshgrid(X,Y);
%N=meshgrid(Y,X);
%plot(X,N,'b');
%plot(M,Y,'b');
%set(gca,'ydir','reverse','xaxislocation','top');
%index = reshape(1:index_X*index_Y,index_Y,index_X);
%for indexi = 1 :index_X
    %for indexj = 1 : index_Y
        %text(length*(indexi-0.5),length*(indexj-0.5),num2str(index(indexj,indexi)));
    %end
%end
set(gca,'YDir','reverse');


