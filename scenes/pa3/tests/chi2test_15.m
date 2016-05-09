obsFrequencies = [ 1, 0, 0, 0, 0, 1, 0, 0, 2, 12, 56, 547, 2829, 2887, 660, 57, 14, 2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 80, 10442, 12026, 117, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13039, 23731, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 545, 22336, 24472, 620, 1, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 225, 4710, 23219, 24078, 5149, 274, 10, 0, 0, 0; 821, 755, 806, 810, 794, 855, 791, 786, 828, 1034, 2776, 10296, 21609, 22028, 10841, 3146, 1031, 820, 842, 758; 2324, 2378, 2415, 2349, 2453, 2359, 2418, 2457, 2716, 3789, 7286, 13930, 20550, 20555, 14416, 7696, 3941, 2807, 2533, 2355; 4055, 4105, 4214, 3969, 4076, 4157, 4347, 4681, 5524, 7536, 11402, 15877, 19280, 19430, 16241, 11444, 7826, 5618, 4557, 4396; 6420, 6249, 6116, 6032, 6105, 6355, 6741, 7641, 9090, 11095, 13997, 17014, 18238, 18612, 16846, 14197, 11417, 9369, 7596, 6901; 10279, 9891, 9659, 9876, 9979, 10260, 10750, 11703, 12701, 13889, 15295, 16244, 17010, 16737, 16192, 15244, 13844, 12839, 11578, 10964 ];
expFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 800.012, 800.002, 800.001, 800.001, 800.002, 800.011, 800.097, 801.253, 818.228, 1031.42, 2833.12, 10322.5, 21440.1, 21779.5, 10818.7, 3009.02, 1057.72, 820.487, 801.41, 800.108; 2404.42, 2401.71, 2401.03, 2401.02, 2401.66, 2404.2, 2415.14, 2468.1, 2729.4, 3855.38, 7396.02, 14132.3, 20320.3, 20481.6, 14463.1, 7633.28, 3946.78, 2752.83, 2472.98, 2416.12; 4115.77, 4065.71, 4048.59, 4048.27, 4064.45, 4112.36, 4241.69, 4592.83, 5508.54, 7599.99, 11344.7, 16008, 19362.3, 19443.4, 16202.5, 11543.5, 7729.68, 5570.87, 4617.76, 4250.93; 6427.13, 6201.89, 6108.14, 6106.25, 6195.43, 6413.37, 6850.05, 7666.72, 9083.93, 11262, 14041.2, 16753.4, 18458.4, 18497.9, 16856.6, 14169.1, 11375, 9163.6, 7714.94, 6876.61; 10335.2, 9986.39, 9819.42, 9815.89, 9975.37, 10315.6, 10868, 11665.4, 12713.1, 13951.7, 15222.5, 16287, 16902.4, 16916.3, 16325.2, 15275.7, 14008.6, 12764.5, 11706.5, 10897.8 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');
