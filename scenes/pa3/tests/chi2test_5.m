obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 48, 2357, 4695, 385, 3, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 63, 7585, 116018, 189243, 32476, 614, 0, 0, 0; 1, 1, 0, 2, 2, 0, 1, 2, 6, 9, 109, 2053, 34896, 213054, 297649, 90472, 7806, 399, 32, 2 ];
expFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 7.99233e-37, 1.60267e-24, 1.9651e-14, 1.85936e-07, 0.000559477, 0.00171477, 1.64615e-05, 5.77738e-11, 8.8427e-20, 3.00598e-31, 0; 8.61362e-36, 0, 0, 0, 0, 0, 0, 1.85559e-39, 6.88864e-31, 3.9478e-22, 4.98056e-14, 2.49266e-07, 0.0138871, 4.17371, 10.0969, 0.324784, 5.43586e-05, 6.82363e-11, 1.83656e-18, 4.92113e-27; 6.23246e-21, 3.99793e-26, 2.61712e-30, 5.85898e-33, 1.74171e-33, 1.09508e-31, 4.15675e-28, 2.73014e-23, 9.08783e-18, 4.27507e-12, 8.36912e-07, 0.0238568, 41.1253, 2370.72, 4688.4, 375.258, 0.91639, 9.90399e-05, 1.04379e-09, 2.82342e-15; 5.69492e-10, 6.21421e-13, 2.81705e-15, 1.04884e-16, 5.62658e-17, 4.97158e-16, 4.73604e-14, 2.54171e-11, 3.78802e-08, 7.68103e-05, 0.110555, 62.0966, 7492.96, 116002, 189469, 32581.7, 631, 2.04531, 0.00198998, 1.05761e-06; 1.72052, 0.85866, 0.546757, 0.427415, 0.413292, 0.488466, 0.677232, 1.22568, 3.03749, 10.3908, 115.103, 2099.43, 34947.4, 212469, 297942, 90542, 7704.21, 389.634, 27.6709, 5.10959 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');
