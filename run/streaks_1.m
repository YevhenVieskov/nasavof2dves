%
%%%%%%%%%%%%%%%%%%%%% MATLAB GRAPHICS %%%%%%%%%%%%%%%%%%%%
%
%  motion picture of streaklines
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  input data
%
icase = input(' Enter case: 1=dam, 2=wave, 3=drop ') 

if icase == 1
  B = load ('dam.str');
  xmin = 0;
  xmax = 10;
  ymin = 0;
  ymax = 5;
elseif icase == 2
  B = load ('wave.str');
  xmin = 0;
  xmax = 40;
  ymin = 0;
  ymax = 10;
else
  B = load ('drop.str');
  xmin = 0;
  xmax = 8;
  ymin = 0;
  ymax = 5;
end
[m,n] = size(B)
%
%  use all the data? or skip some values?
%
skip = input(' Enter SKIP to use every SKIP-th value of the data ') 
%
% determine the number of frames and the size of each
%
frames = 0;
for i=1:1:m
   if B(i,2) == -1 & B(i,3) == -1 & B(i,4) == -1
      frames = frames + 1
      count(frames) = 0;
   else
      count(frames) = count(frames) + 1; 
   end
end
%
still = input(' Enter the desired frame number ')
%
index = 0;
for j = 1:1:frames
   index = index + 1;
   for i=1:1:count(j)/skip
      index = index + skip;
      if (j == still)
         x(i) = B(index,1);
         y(i) = B(index,2);
      end
   end
end
plot(x,y,'.');
set(gca,'Xlim',[xmin xmax]);
set(gca,'Ylim',[ymin ymax]);
xlabel('x','FontName','Helvetica','FontWeight','bold','FontSize',16);
ylabel('y','FontName','Helvetica','FontWeight','bold','FontSize',16,'Rotation',0);
title('Streaklines','FontName','Helvetica','FontWeight','bold','FontSize',16);
