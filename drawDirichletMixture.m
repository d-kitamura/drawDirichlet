% Drawing three-dimensional Dirichlet distribution 
% using standard simplex triangle

clear; close all; clc;

%% Set parameters
nBin = 200; % number of bins for discretization
alpha = [ % parameter of Dirichlet distribution mixture (column vector corresponds to parameter vector)
    0.5, 1, 2, 3;
    0.5, 2, 4, 1;
    0.5, 1, 8, 3;
    ];
gamma = [0.3, 0.4, 0.2, 0.1]; % contribution ratioc

%% Define axis
dim = 3; % dimension (only dim=3 is supported)
nDist = size(alpha, 2); % number of dirichlet distributions
axisX = linspace(0, 1, nBin); % discrete axis
axisY = linspace(0, 1, nBin); % discrete axis
[ax, ay] = meshgrid(axisX, axisY); % two-dimensional discrete space
ax = ax .* fliplr(tril(ones(length(ax)), 1)');
ay = ay .* fliplr(tril(ones(length(ay)), 1)');
z = 1 - ax - ay;
z = z .* fliplr(tril(ones(length(z)), -1)');

%% Calculate PDF value
v = zeros(nBin, nBin, nDist);
for i = 1:nDist
    v(:, :, i) = alpha(i) * calcDirishletPdfVal(dim, alpha(:, i), ax, ay, z);
end
pdfVal = sum(v, 3);

%% Affine transform
M = [1, 1/2;
    0, sqrt(3)/2];
tmp = cat(1, permute(ax, [3, 1, 2]), permute(ay, [3, 1, 2]));
axisSimplex = pagemtimes(M, tmp);
axSimplex = squeeze(axisSimplex(1, :, :));
aySimplex = squeeze(axisSimplex(2, :, :));

%% Show 3D graph
zMaxCond = 15;
fig1 = figure("Position", [50, 50, 561, 473]);
surf(axSimplex, aySimplex, pdfVal, "EdgeColor", "none");
xlim([0, 1]);
ylim([0, sqrt(3)/2]);
maxVal = max(pdfVal(pdfVal~=Inf), [], "all");
if maxVal >= zMaxCond % 最大でもzMaxCondとし，それ以外は最大値の1.2倍
    zMax = zMaxCond; % ここは要調整
    cMax = zMaxCond; % カラーバーの上限値
else
    zMax = ceil(maxVal * 1.2); % ここは要調整
    cMax = maxVal; % カラーバーの上限値
end
zlim([0, zMax]); clim([0, cMax]);
grid on; box on;
view([-15, 45]); % 3D回転の視点
xticks([]); yticks([]); % 軸のメモリを削除
zlabel("Probabilistic density", ...
    "FontSize", 12, ...
    "FontName", "Arial");
set(gca, "FontSize", 12); % 軸の値のフォントサイズ
str = "$$\alpha = [" + alpha(1) + ", " + alpha(2) + ", " + alpha(3) + "]$$";
text(0.9, 0.9, 12, str, ...
    "FontSize", 14, ...
    "EdgeColor", [0, 0, 0], ...
    "BackgroundColor", [1, 1, 1], ...
    "Interpreter", "latex");

%% Show top-view graph
fig2 = figure("Position", [50, 50, 561, 473]);
surf(axSimplex, aySimplex, pdfVal, "EdgeColor", "none");
xlim([0, 1]);
ylim([0, sqrt(3)/2]);
maxVal = max(pdfVal(pdfVal~=Inf), [], "all");
zMax = maxVal;
cMax = maxVal/4; % カラーバーの上限値（ここは要調整）
zlim([0, zMax]); clim([0, cMax]);
grid on; box on;
view([0, 90]); % 3D回転の視点
xticks([]); yticks([]); % 軸のメモリを削除
zlabel("Probabilistic density", ...
    "FontSize", 12, ...
    "FontName", "Arial");
set(gca, "FontSize", 12); % 軸の値のフォントサイズ
str = "$$\alpha = [" + alpha(1) + ", " + alpha(2) + ", " + alpha(3) + "]$$";
text(0.65, 0.8, 0, str, ...
    "FontSize", 14, ...
    "EdgeColor", [0, 0, 0], ...
    "BackgroundColor", [1, 1, 1], ...
    "Interpreter", "latex");

%% Save figures
savefig(fig1, "./simplex_3d.fig");
saveas(fig1, "./simplex_3d.svg", "svg");
saveas(fig1, "./simplex_3d.png", "png");

savefig(fig2, "./simplex_topview.fig");
saveas(fig2, "./simplex_topview.svg", "svg");
saveas(fig2, "./simplex_topview.png", "png");

%% Local function
function v = calcDirishletPdfVal(dim, alpha, ax, ay, z)
% g = @(alpha, d) (gamma(sum(alpha)) / (prod(gamma(alpha)))); % 多変量ベータ関数
g = @(alpha, d) exp(gammaln(sum(alpha)) - sum(gammaln(alpha))); % 多変量ベータ関数（オーバーフロー対策）
v = g(alpha, dim) .* ax .^ (alpha(1) - 1) .* ay .^ (alpha(2) - 1) .* z .^ (alpha(3) - 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%