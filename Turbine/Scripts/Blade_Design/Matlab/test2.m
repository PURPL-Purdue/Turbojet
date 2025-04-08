clf

t = 0:8;
xy = [0 0;1 1; 1.7 0;1 -1;0 0; -1 1; -1.7 0; -1 -1; 0 0].';
infty = csape(t, xy, 'periodic');
diff = fnder(infty,2);

hold on
fnplt(diff, 2)
fnplt(infty, 2)
axis([-2 2 -1.1 1.1])
plot(xy(1,:),xy(2,:),'o')
hold off