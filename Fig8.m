clear;
close all;

ii=0;
jj=0;
 for al=0:0.01:1
     jj=0;
     ii=ii+1;
     for th=0:0.01:pi/2
         jj=jj+1;
         fid(ii,jj)=(al^2*cos(2*th) - cos(2*th)/2 + al*sin(2*th)*(1 - al^2)^(1/2) + 1/2)/(2*cos(2*th)*al^2 - cos(2*th) + 1);
     end
 end
 surf(fid)
 %y=0:0.01:1;
 %x=0:0.01:pi/2;
% contour(x,y,fid)
 ylabel('\alpha')
 xlabel('\theta')
 xticks([0, 39.25, 78.5, 117.75, 157])
 xticklabels({0, '\pi/8', '\pi/4', '3\pi/8', '\pi/2'})
 yticks([0,20,40,60,80])
 yticklabels({0,'0.2','0.4','0.6','0.8'})
 set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 12)
set(gca, 'FontWeight', 'bold')
 axis tight
 view(0, 90)
 colorbar