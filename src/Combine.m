data1 = load("../Data/HTBData/HTB2-mCitrine_YPD_M2.mat");
data2 = load("../Data/HTBData/HTB2_YPD_parameters_full.mat");

pulsedata = struct;

a = data1.pulsedata.totvolG1_daughters;
b = data2.pulsedata.totvolG1_daughters;
sz1 = size(a);
sz2 = size(b);
dim = min(sz1(2), sz2(2));
pulsedata.totvolG1_daughters = cat(1, a(:, 1:dim), b(:, 1:dim));

a = data1.pulsedata.volumestartG1_daughters;
b = data2.pulsedata.volumestartG1_daughters;
pulsedata.volumestartG1_daughters = cat(2, a, b);

a = data1.pulsedata.volumeendG1_daughters;
b = data2.pulsedata.volumeendG1_daughters;
pulsedata.volumeendG1_daughters = cat(2, a, b);

a = data1.pulsedata.volumeG1_daughters;
b = data2.pulsedata.volumeG1_daughters;
sz1 = size(a);
sz2 = size(b);
dim = min(sz1(2), sz2(2));
pulsedata.volumeG1_daughters = cat(1, a(:, 1:dim), b(:, 1:dim));

a = data1.pulsedata.totvolSG2M;
b = data2.pulsedata.totvolSG2M;
sz1 = size(a);
sz2 = size(b);
dim = min(sz1(2), sz2(2));
pulsedata.totvolSG2M = cat(1, a(:, 1:dim), b(:, 1:dim));

a = data1.pulsedata.volbud;
b = data2.pulsedata.volbud;
sz1 = size(a);
sz2 = size(b);
dim = min(sz1(2), sz2(2));
pulsedata.volbud = cat(1, a(:, 1:dim), b(:, 1:dim));

save("CombinedYPD", "pulsedata");