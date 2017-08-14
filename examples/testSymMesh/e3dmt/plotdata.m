
d1 = load('fields_1_1.txt');
d2 = load('fields_1_2.txt');

x = d1(:,1);
y = d1(:,2);


d1(:,4:9) = -d1(:,4:9);
d2(:,4:9) = -d2(:,4:9);

subplot(2,6,1)
d = reshape(d1(:,4), 3,3);
imagesc(d)
title("Hx r")
colorbar

subplot(2,6,2)
d = reshape(d1(:,5), 3,3);
imagesc(d)
title("Hx i")
colorbar

subplot(2,6,3)
d = reshape(d1(:,6), 3,3);
imagesc(d)
title("Hy r")
colorbar

subplot(2,6,4)
d = reshape(d1(:,7), 3,3);
imagesc(d)
title("Hy i")
colorbar

subplot(2,6,5)
d = reshape(d1(:,8), 3,3);
imagesc(d)
title("Hz r")
colorbar

subplot(2,6,6)
d = reshape(d1(:,9), 3,3);
imagesc(d)
title("Hz i")
colorbar

%-----------

subplot(2,6,7)
d = reshape(d2(:,4), 3,3);
imagesc(-d')
title("Hx r")
colorbar

subplot(2,6,8)
d = reshape(d2(:,5), 3,3);
imagesc(-d')
title("Hx i")
colorbar

subplot(2,6,9)
d = reshape(d2(:,6), 3,3);
imagesc(-d')
title("Hy r")
colorbar

subplot(2,6,10)
d = reshape(d2(:,7), 3,3);
imagesc(-d')
title("Hy i")
colorbar

subplot(2,6,11)
d = reshape(d2(:,8), 3,3);
imagesc(-d')
title("Hz r")
colorbar

subplot(2,6,12)
d = reshape(d2(:,9), 3,3);
imagesc(-d')
title("Hz i")
colorbar
