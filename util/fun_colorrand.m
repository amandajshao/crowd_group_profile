function crand = fun_colorrand(index)

%% construct enough colors for different groups
% colornum = 100;
% cmap = hsv(colornum);
% rr = rgb2hsv([1,0,0]);
% colormatrix = cmap(randperm(length(cmap(:,1))), :);
% load colormatrix3.mat
% crand = colormatrix(index,:);

%%
c1 = [255, 128, 255];
c2 = [0, 128, 255];
c3 = [0, 255, 128];
c4 = [255, 128, 64];
c5 = [128, 0, 255];
c6 = [0, 128, 0];
c7 = [255, 0, 128];
c8 = [255, 255, 128];
c9 = [128, 255, 255];
c10 = [100, 0, 100]; 
c11 = [0, 128, 128];
c12 = [128, 128, 192];
c13 = [0, 0, 215];
c14 = [255, 70, 20];
c15 = [128, 128, 0];
c16 = [40, 40, 200];
c17 = [20, 180, 30];
c18 = [230, 200, 16];
c19 = [220, 0 ,220];
c20 = [255, 200, 150];
c21 = [155, 220, 90];
c22 = [55, 155, 255];
c23 = [183,0,0];
c24 = [255,255,0];
c25 = [240,120,140];
c26 = [240,70,10];
c27 = [48,203,196];
c28 = [128, 128, 128];
c29 = [0,128,128];
c30 = [102, 0, 51];


colormatrix = [c1;c2;c3;c4;c5;c6;c7;c8;c9;c10;...
    c11;c12;c13;c14;c15;c16;c17;c18;c19;c20;...
    c21;c22;c23;c24;c25;c26;c27;c28;c29;c30];
colormatrix = colormatrix./255;
crand = colormatrix(index,:);

end