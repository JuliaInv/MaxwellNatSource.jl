
d1 = load("data1.txt");
d2 = load("data2.txt");

d3 = load("data.csv");

d2 = d2(2:end-1, :);

% 500 Hz
%xI = d3(:,27) ;  xQ = d3(:,28);
%zI = d3(:,107) ; zQ = d3(:,108);

% 12 Hz
%xI = d3(:,7) ;  xQ = d3(:,8);
%zI = d3(:,87) ; zQ = d3(:,88);

% 19905 Hz
xI = d3(:,9) ;  xQ = d3(:,10);
zI = d3(:,89) ; zQ = d3(:,90);

figure(1)

subplot(2,2, 1)
hold off
plot(d3(:,1), xI, 'r')
hold on
plot(d3(:,1), xQ, 'b')
axis tight
legend('SCX19905I','SCX19905Q')


subplot(2,2, 2)
hold off
plot(d3(:,1), zI, 'r')
hold on
plot(d3(:,1), zQ, 'b')
axis tight
legend('SCZ19905I','SCZ19905Q')



subplot(2,2, 3)
hold off
plot(d2(:,1), d2(:,4), 'r')
hold on
plot(d2(:,1), d2(:,5)-3.41e-5, 'b')
legend('r','i-3.41e-5')
axis tight
title('polarization 2, Hx')


subplot(2,2, 4)
hold off
plot(d2(:,1), d2(:,8), 'r')
hold on
plot(d2(:,1), d2(:,9), 'b')
legend('r','i')
axis tight
title('polarization 2, Hz')
