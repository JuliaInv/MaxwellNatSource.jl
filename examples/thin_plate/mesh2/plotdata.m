

d3 = load("../data.csv");

% 500 Hz
%xI = d3(:,27) ;  xQ = d3(:,28);
%yI = d3(:,67) ;  yQ = d3(:,68);
%zI = d3(:,107) ; zQ = d3(:,108);

% 12 Hz
xI = d3(:,7) ;  xQ = d3(:,8);
yI = d3(:,47) ;  yQ = d3(:,48);
zI = d3(:,87) ; zQ = d3(:,88);


d1 = load("fields_2_1.txt");
d2 = load("fields_2_2.txt");

d1hs = load("halfspace/fields_2_1.txt");
d2hs = load("halfspace/fields_2_2.txt");


d2 = d2(2:end-1, :);
d1 = d1(2:end-1, :);

d2hs = d2hs(2:end-1, :);
d1hs = d1hs(2:end-1, :);


d1(:,4:end) = -(d1(:,4:end) - d1hs(:,4:end));
d2(:,4:end) = -(d2(:,4:end) - d2hs(:,4:end));


d1(:,4:7) = -d1(:,4:7);
d2(:,4:7) = -d2(:,4:7);

figure(1)



subplot(3,3, 1)
hold off
plot(d1(:,1), d1(:,4), 'r')
hold on
plot(d1(:,1), d1(:,5), 'b')
axis tight
legend('r','i')
title("Hx")


subplot(3,3, 2)
hold off
plot(d1(:,1), d1(:,6), 'r')
hold on
plot(d1(:,1), d1(:,7), 'b')
axis tight
legend('r','i')
title("Hy")

subplot(3,3, 3)
hold off
plot(d1(:,1), d1(:,8), 'r')
hold on
plot(d1(:,1), d1(:,9), 'b')
axis tight
legend('r','i')
title("Hz")


%---

subplot(3,3, 4)
hold off
plot(d1(:,1), d2(:,4), 'r')
hold on
plot(d1(:,1), d2(:,5), 'b')
axis tight
legend('r','i')
title("Hx")


subplot(3,3, 5)
hold off
plot(d1(:,1), d2(:,6), 'r')
hold on
plot(d1(:,1), d2(:,7), 'b')
axis tight
legend('r','i')
title("Hy")

subplot(3,3, 6)
hold off
plot(d1(:,1), d2(:,8), 'r')
hold on
plot(d1(:,1), d2(:,9), 'b')
axis tight
legend('r','i')
title("Hz")

%---

subplot(3,3, 7)
hold off
plot(d3(:,1), xI, 'r')
hold on
plot(d3(:,1), xQ, 'b')
axis tight
legend('SCXI','SCXQ')

subplot(3,3, 8)
hold off
plot(d3(:,1), yI, 'r')
hold on
plot(d3(:,1), yQ, 'b')
axis tight
legend('SCYI','SCYQ')

subplot(3,3, 9)
hold off
plot(d3(:,1), zI, 'r')
hold on
plot(d3(:,1), zQ, 'b')
axis tight
legend('SCZI','SCZQ')


