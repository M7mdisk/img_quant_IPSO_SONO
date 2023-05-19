clear all
close all
clc

pkg load statistics

ps = 10;
nfe_max = 100;
Xmin = 0;
Xmax = 255;
K = 12;
D = K * 3;
fhd=@mymse;

file_name = "Mona_Lisa";
ext = ".jpg";

originalImage = imread(strcat(file_name,ext));

[gbestval,ccurve, dcurve,gbest] = PSO_sono_CEC2022(ps, nfe_max, Xmin, Xmax, D,fhd, originalImage);
[gbestvalx,ccurvex, dcurvex,gbestx] = IPSO_sono_CEC2022(ps, nfe_max, Xmin, Xmax, D,fhd, originalImage);

disp(["Regular PSO fit:" num2str(gbestval)]);
disp(["IPSO fit:",num2str(gbestvalx)]);
paletteToImg(gbest,originalImage,sprintf("output/%s_%d_%d.jpg",file_name,K,nfe_max));
paletteToImg(gbestx,originalImage,sprintf("output/%s_%d_%d_IMP.jpg",file_name,K,nfe_max));


