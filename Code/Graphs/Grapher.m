clear;
kratos=readtable('KratosCylinderVisc0.1.txt');
kratos=table2array(kratos);
plot(kratos(:,1),kratos(:,2))
mine2=readtable("MineCylinderVisc0.0816.txt");
mine2=table2array(mine2);
hold on
plot(mine2(:,1),mine2(:,2))
legend('Kratos','Mine')
figure
kratos=readtable('LidDrivenNSVisc0.1Kratos.txt');
kratos=table2array(kratos);
plot(kratos(:,1),kratos(:,2))
mine=readtable("LidDriven75NSVisc0.081633Mine.txt");
mine=table2array(mine);
hold on
plot(mine(:,1),mine(:,2))
mine2=readtable("LidDriven20NSVisc0.081633Mine.txt");
mine2=table2array(mine2);
hold on
plot(mine2(:,1),mine2(:,2))
legend('Kratos','Mine75','Mine20')