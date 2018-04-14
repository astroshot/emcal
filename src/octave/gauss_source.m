% % generate gauss beam
% % freq,ds,w0,z
% % 生成gauss波束
% % freq为频率/GHz
% % r为网格距离矩阵
% % w0为束腰半径m
% % z为待观察平面
function source = gauss_source(freq, w0, ds, N, Z)

    C0 = 299792458;  %光速
    wavelength = C0 / freq / 1e9;
    k = 2 * pi / wavelength;

    x = ((1:N) - N / 2) * ds;
    y = x;
    [Y, X] = meshgrid(x, y);
    r = sqrt(X .^ 2 + Y .^ 2 + Z .^ 2);

    w = w0 * sqrt(1 + (wavelength * Z / pi / w0 / w0) ^ 2);
    Rdl = Z / (Z .^ 2 + (pi * w0 * w0 / wavelength) .^ 2);
    phi = atan(wavelength * Z / pi / w0 / w0);
    complex = -r .^ 2 ./ (w ^ 2) - 1j * k * Z - 1j * pi * r .^ 2 * Rdl / wavelength + 1j .* phi;
    source = sqrt(2 / pi ./(w .^ 2)) .* exp(complex);

end
