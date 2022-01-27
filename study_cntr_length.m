ccc



g = ones(100);

g(12:27,23:48) = 2;
g(15:20,30:45) = 3;

x = 1:100;
y = 101:200;

cntr_plot = [1:3];
cntr = [2];

[Cf,ch] = contourf(x,y,g,cntr_plot);

hold on
[C,h] = contour(x,y,...
    g,[1 1]*cntr,'w','linewidth',2);
C(:,C(1,:)<11) = [];

d = zeros(length(cntr),1);
Vq = interp2(x,y,g,C(1,:),C(2,:));
j  =1;
for i = 1:length(cntr)
    %         bM = regionprops(meanSSH>=cntr(i),'Perimeter');
    %         cntr_length(i,j) = bM.Perimeter;
    inds = find(abs(Vq-cntr(i))<1e-10);
    for k = 1:length(inds)-1
        delta = norm(C(:,inds(k+1))-C(:,inds(k)));
        if delta < 2
            d(i,j) = d(i,j) + delta;
        end
    end
end









