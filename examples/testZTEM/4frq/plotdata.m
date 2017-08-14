

d1 = load("dataZTEM_8.txt");
d2 = load("dataZTEM_31.txt");
d3 = load("dataZTEM_199.txt");
d4 = load("dataZTEM_1255.txt");

%d_e3d = load("e3dmt/MT_data.txt");
d_e3d = load("e3dmt/dpred0.txt");


nfrq = 4;
d_e3d = reshape(d_e3d, 625, nfrq, 7);


d = d4;
d_e3d = d_e3d(:,4,:);



x = d(:,1);
y = d(:,2);

S = 30;
figure(1)

subplot(2,2, 1)
scatter(x,y, S, d(:,4), 'filled')
axis image
title("Tx real")
colorbar

subplot(2,2, 2)
scatter(x,y, S, d(:,5), 'filled')
axis image
title("Tx imaginary")
colorbar

subplot(2,2, 3)
scatter(x,y, S, d(:,6), 'filled')
axis image
title("Ty real")
colorbar

subplot(2,2, 4)
scatter(x,y, S, d(:,7), 'filled')
axis image
title("Ty imaginary")
colorbar

%---------------------------------------------

figure(2)

subplot(2,2, 1)
scatter(x,y, S, d_e3d(:,4), 'filled')
axis image
title("Tx real")
colorbar

subplot(2,2, 2)
scatter(x,y, S, d_e3d(:,5), 'filled')
axis image
title("Tx imaginary")
colorbar

subplot(2,2, 3)
scatter(x,y, S, d_e3d(:,6), 'filled')
axis image
title("Ty real")
colorbar

subplot(2,2, 4)
scatter(x,y, S, d_e3d(:,7), 'filled')
axis image
title("Ty imaginary")
colorbar
