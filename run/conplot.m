%
%%%%%%%%%%%%%%%% MATLAB CONTOUR GRAPHICS %%%%%%%%%%%%%%%%%
%
% contour plot template for discrete data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
B = load ('contour_data');
[m,n] = size(B)

nx = input(' Enter the number of grid intervals in x ');
ny = input(' Enter the number of grid intervals in y ');

index = 0;
for j=1:1:ny+1      % no. gridpoints in y
   for i=1:1:nx+1   % no. gridpoints in x
      index = index + 1;
      X(i,j) = B(index,1);
      Y(i,j) = B(index,2);
      Z(i,j) = B(index,3);
   end
end

no_levels = 15;
[C1,h1] = contourf(X,Y,Z,no_levels);
colormap(jet);

xlabel('x','FontName','Helvetica','FontWeight','bold','FontSize',16);
ylabel('z','FontName','Helvetica','FontWeight','bold','FontSize',16,'Rotation',0);
title('streamfunction','FontName','Helvetica','FontWeight','bold','FontSize',16);
