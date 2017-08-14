
d = load('dpred0.txt');

nd = size(d,1);

sd = zeros(nd,4);
sd = abs(d(:,4:7))*0.05 + 1.e-4;


datainv = zeros(nd,3+4*2);

datainv(:,1:3) = d(:,1:3);
datainv(:,4:2:end) = d(:,4:7);
datainv(:,5:2:end) = sd;

save data_inv.txt datainv -ascii
