obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 50565, 98, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 299, 390848, 4962, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2940, 466503, 18968, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1339, 57302, 5104, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 101, 24, 0, 0, 0, 0, 0, 0, 0 ];
expFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.38436e-21, 1.56638, 50124.3, 76.2399, 1.74459e-16, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 3.72011e-30, 8.05962e-10, 301.798, 391146, 5028.03, 4.5545e-07, 6.69325e-26, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1.03057e-31, 3.39362e-16, 0.000159113, 2891.51, 466974, 18833.7, 0.00721161, 8.62208e-14, 1.09617e-28, 0, 0, 0, 0; 0, 0, 0, 0, 0, 4.98355e-37, 3.60358e-27, 1.04754e-17, 1.94658e-09, 0.0145155, 1385.92, 57085.1, 5049.47, 0.171486, 4.29408e-08, 4.08291e-16, 2.07047e-25, 2.86645e-35, 0, 0; 1.56503e-13, 1.56503e-13, 1.56503e-13, 1.56503e-13, 1.56503e-13, 1.56503e-13, 1.57633e-13, 9.9474e-12, 9.15259e-08, 0.00537163, 9.78879, 118.624, 21.9659, 0.0292074, 6.83308e-07, 4.65177e-11, 1.62543e-13, 1.56503e-13, 1.56503e-13, 1.56503e-13 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');
