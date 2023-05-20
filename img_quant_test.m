clear all
close all
clc

pkg load statistics

Xmin = 0;
Xmax = 255;

% Population size
ps = 100;

% Max function evals
nfe_max = 1000;

% Number of colors to be used

fhd=@mymse;

file_name = "Mona_Lisa";
ext = ".jpg";

originalImage = imread(strcat(file_name,ext));
if size(originalImage, 3) == 4
    originalImage = originalImage(:, :, 1:3);
end

for K = [5 7 10 12 15 20]
  D = K * 3;
  [gbestval,ccurve, dcurve,gbest] = PSO_sono_CEC2022(ps, nfe_max, Xmin, Xmax, D,fhd, originalImage);
  [gbestvalx,ccurvex, dcurvex,gbestx] = IPSO_sono_CEC2022(ps, nfe_max, Xmin, Xmax, D,fhd, originalImage);
  disp(["Number of Colors:",num2str(K)]);
  disp(["Regular PSO fit:" num2str(gbestval)]);
  disp(["IPSO fit:",num2str(gbestvalx)]);
  paletteToImg(gbest,originalImage,sprintf("output/%s_%d_%d.jpg",file_name,K,nfe_max));
  paletteToImg(gbestx,originalImage,sprintf("output/%s_%d_%d_IMP.jpg",file_name,K,nfe_max));
endfor

