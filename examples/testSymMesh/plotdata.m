
d1 = load('data11.txt');
d2 = load('data22.txt');

dZ = load('dataZTEM.txt');

x = d1(:,1);
y = d1(:,2);


subplot(3,6,1)
d = reshape(d1(:,4), 7,7);
imagesc(d)
title("Hx r")
colorbar

subplot(3,6,2)
d = reshape(d1(:,5), 7,7);
imagesc(d)
title("Hx i")
colorbar

subplot(3,6,3)
d = reshape(d1(:,6), 7,7);
imagesc(d)
title("Hy r")
colorbar

subplot(3,6,4)
d = reshape(d1(:,7), 7,7);
imagesc(d)
title("Hy i")
colorbar

subplot(3,6,5)
d = reshape(d1(:,8), 7,7);
imagesc(d)
title("Hz r")
colorbar

subplot(3,6,6)
d = reshape(d1(:,9), 7,7);
imagesc(d)
title("Hz i")
colorbar

%-----------

subplot(3,6,7)
d = reshape(d2(:,4), 7,7);
imagesc(d')
title("Hx r")
colorbar

subplot(3,6,8)
d = reshape(d2(:,5), 7,7);
imagesc(d')
title("Hx i")
colorbar

subplot(3,6,9)
d = reshape(d2(:,6), 7,7);
imagesc(d')
title("Hy r")
colorbar

subplot(3,6,10)
d = reshape(d2(:,7), 7,7);
imagesc(d')
title("Hy i")
colorbar

subplot(3,6,11)
d = reshape(d2(:,8), 7,7);
imagesc(d')
title("Hz r")
colorbar

subplot(3,6,12)
d = reshape(d2(:,9), 7,7);
imagesc(d')
title("Hz i")
colorbar

%-----------

subplot(3,6,13)
d = reshape(dZ(:,4), 7,7);
imagesc(d)
title("T1 r")
colorbar

subplot(3,6,14)
d = reshape(dZ(:,5), 7,7);
imagesc(d)
title("T1 i")
colorbar

subplot(3,6,15)
d = reshape(dZ(:,6), 7,7);
imagesc(d)
title("T2 r")
colorbar

subplot(3,6,16)
d = reshape(dZ(:,7), 7,7);
imagesc(d)
title("T2 i")
colorbar

