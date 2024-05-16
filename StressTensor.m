clear
clc
close

format longG

components = {-22.2 9.1 7.3 9.1 -16.9 -4.6 7.3 -4.6 31.8};
[sigma_x, tau_xy, tau_xz, tau_yx, sigma_y, tau_yz, tau_zx, tau_zy, sigma_z] = components{:};

% symbolic
% components = {-22.2 9.1 7.3 9.1 -16.9 -4.6 7.3 -4.6 31.8};
% components_sym = cellfun(@sym, components, 'UniformOutput', 0);
% [sigma_x, tau_xy, tau_xz, tau_yx, sigma_y, tau_yz, tau_zx, tau_zy, sigma_z] = components_sym{:};

bf_sigma = [sigma_x tau_xy tau_xz;
           tau_yx sigma_y tau_yz;
           tau_zx tau_zy sigma_z]

P = -(sigma_x/3 + sigma_y/3 + sigma_z/3)
bf_sigma_hydrostatic = -P*eye(3)
bf_sigma_prime = bf_sigma - bf_sigma_hydrostatic

I_1 = sigma_x + sigma_y + sigma_z
I_2 = sigma_x*sigma_y + sigma_x*sigma_z + sigma_y*sigma_z - tau_xy*tau_yx - tau_xz*tau_zx - tau_yz*tau_zy
I_3 = sigma_x*sigma_y*sigma_z - sigma_z*tau_xy*tau_yx - sigma_y*tau_xz*tau_zx - sigma_x*tau_yz*tau_zy + tau_xy*tau_zx*tau_yz + tau_yx*tau_xz*tau_zy

% % via bf_sigma
% J_1 = 0
% J_2 = sigma_x^2/3 - (sigma_x*sigma_y)/3 - (sigma_x*sigma_z)/3 + sigma_y^2/3 - (sigma_y*sigma_z)/3 + sigma_z^2/3 + tau_xy*tau_yx + tau_xz*tau_zx + tau_yz*tau_zy
% J_3 (2*sigma_x^3)/27 - (sigma_x^2*sigma_y)/9 - (sigma_x*sigma_z^2)/9 - (sigma_x^2*sigma_z)/9 - (sigma_y*sigma_z^2)/9 - (sigma_y^2*sigma_z)/9 - (sigma_x*sigma_y^2)/9 + (2*sigma_y^3)/27 + (2*sigma_z^3)/27 + (4*sigma_x*sigma_y*sigma_z)/9 + (sigma_x*tau_xy*tau_yx)/3 + (sigma_y*tau_xy*tau_yx)/3 + (sigma_x*tau_xz*tau_zx)/3 - (2*sigma_z*tau_xy*tau_yx)/3 - (2*sigma_y*tau_xz*tau_zx)/3 - (2*sigma_x*tau_yz*tau_zy)/3 + (sigma_z*tau_xz*tau_zx)/3 + (sigma_y*tau_yz*tau_zy)/3 + (sigma_z*tau_yz*tau_zy)/3 + tau_xy*tau_zx*tau_yz + tau_yx*tau_xz*tau_zy;

J_1 = 0
J_2 = 1/3*I_1^2 - I_2
J_3 = 2/27*I_1^3 - 1/3*I_1*I_2 + I_3

theta = asin(-J_3/2*(3/J_2)^(3/2))/3

sigma_1 = I_1/3 + sqrt(J_2)*(cos(theta) - 1/sqrt(3)*sin(theta))
sigma_2 = I_1/3 + 2/sqrt(3)*sqrt(J_2)*sin(theta)
sigma_3 = I_1/3 - sqrt(J_2)*(cos(theta) + 1/sqrt(3)*sin(theta))

% sorted_pr = cellfun(@(x) sort(x, "descend"), {sigma_1, sigma_2, sigma_3}, "UniformOutput", 0)
% [sigma_1, sigma_2, sigma_3] = sorted_pr{:}

sigma_1_prime = sigma_1 - I_1/3 
sigma_2_prime = sigma_2 - I_1/3 
sigma_3_prime = sigma_3 - I_1/3

tau_13 = 1/2*(sigma_1 - sigma_3)
tau_12 = 1/2*(sigma_1 - sigma_2)
tau_23 = 1/2*(sigma_2 - sigma_3)

sigma_13 = 1/2*(sigma_1 + sigma_3)
sigma_12 = 1/2*(sigma_1 + sigma_2)
sigma_23 = 1/2*(sigma_2 + sigma_3)

sigma_oct = 1/3*I_1

% tau_oct = 1/3*sqrt((sigma_1 - sigma_2)^2 + (sigma_2 - sigma_3)^2 + (sigma_3 - sigma_1)^2)
tau_oct = sqrt(2/3*J_2)

bar_bf_sigma_prime = bf_sigma_prime/tau_oct

% sigma_vM = sqrt(sigma_x^2 - sigma_x*sigma_y - sigma_x*sigma_z + sigma_y^2 - sigma_y*sigma_z + sigma_z^2 + 3*tau_xy*tau_yx + 3*tau_xz*tau_zx + 3*tau_yz*tau_zy)
% sigma_vM = sqrt(sigma_1^2 - sigma_1*sigma_2 - sigma_1*sigma_3 + sigma_2^2 - sigma_2*sigma_3 + sigma_3^2)
sigma_vM = sqrt(3*J_2)

figure()
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

grid on;
box on;
hold on;
axis equal;
axis tight;
xlabel('$\sigma_n$', 'FontSize', 14);
ylabel('$\tau_n$', 'FontSize', 14);
title("The Mohr's Diagram", 'FontSize', 16);

double_pr = cellfun(@double, {sigma_1, sigma_2, sigma_3}, "UniformOutput", 0)
[sigma_1_double, sigma_2_double, sigma_3_double] = double_pr{:}

R_1 = tau_23
R_2 = tau_13
R_3 = tau_12

O_1 = [sigma_23, 0]
O_2 = [sigma_13, 0]
O_3 = [sigma_12, 0]

angle = (0:pi/1800:pi);

x_1 = O_1(1) * ones([1 1801]) + R_1 * cos(angle);
y_1 = R_1 * sin(angle);

x_2 = O_2(1) * ones([1 1801]) + R_2 * cos(angle);
y_2 = R_2 * sin(angle);

x_3 = O_3(1) * ones([1 1801]) + R_3 * cos(angle);
y_3 = R_3* sin(angle);

plot(x_1, y_1, '-k', "LineWidth", 1.5)
plot(x_2, y_2, '-k', "LineWidth", 1.5)
plot(x_3, y_3, '-k', "LineWidth", 1.5)
fill([x_1, x_2, x_3], [y_1, y_2, y_3], 'k', 'FaceAlpha', 0.5);

% exportgraphics(gcf, "The_Mohrs_Diagram.pdf", "ContentType","vector")
% exportgraphics(gcf, "The_Mohrs_Diagram.eps", "ContentType","vector")
% exportgraphics(gcf, "The_Mohrs_Diagram.png", "Resolution", 1200)
