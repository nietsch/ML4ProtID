from pyopenms import *
#import os
#import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#tol = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0]
#tol = [float(t) for t in range(5, 151)]
#tol = np.logspace(1.0, 3.0, num=100)
#tol = np.logspace(1.0, 4.0, num=50)
#tol = [13.0]

filenames = ["1_A", "1_B", "1_C", "2_A", "2_B", "2_C", "3_A", "3_B", "3_C", "4_A", "4_B", "4_C"]

res_only_filtering = [18709,19433,17434,17435,18693,17571,17580,18880,18388,16754,18049,17048]
res_only_perc = [19085,19891,17838,17778,19154,17942,17961,19185,18700,17127,18359,17449]
res_SA = [19286,19982,17940,17865,19238,18200,18129,19402,18783,17240,18568,17622]
res_RTdiff = [19701,20483,18577,18577,19673,18845,18571,19686,19291,17870,19114,18163]
res_SA_RTdiff = [19908,20697,18845,18650,19902,18937,18915,19904,19513,17993,19336,18426]
res_SA_RTdiff_absRT = [20006,20721,18858,18711,19908,19062,18966,19982,19624,18038,19365,18441]
res_RTdiff_absRT = [19721,20579,18657,18607,19680,18884,18778,19713,19363,17946,19203,18255]


res_only_filtering_top10 = [18709,19433,17434,17435,18693,17571,17580,18880,18388,16754,18049,17048]
res_only_perc_top10 = [18442,19168,17206,17153,18574,17271,17409,18568,18100,16463,17735,16795]
res_SA_top10 = [18079,18772,16967,16784,17935,17348,17171,18394,17726,16054,17807,16708]
res_RTdiff_top10 = [18777,19335,18034,17546,18891,17654,17794,19033,18449,16863,18476,17238]
res_SA_RTdiff_top10 = [18407,19009,17341,17044,18509,17187,17323,18549,18193,16687,18032,17150]
res_SA_RTdiff_absRT_top10 = [18430,19093,17601,17170,18441,17409,17367,18530,18190,16716,18056,16918]
res_RTdiff_absRT_top10 = [18852,19769,17751,17895,18985,18234,17969,19069,18675,16782,18526,17483]


lst1 = [res_only_filtering[0], res_only_perc[0], res_SA[0], res_RTdiff[0], res_SA_RTdiff[0], res_SA_RTdiff_absRT[0]]
lst2 = [res_only_filtering[1], res_only_perc[1], res_SA[1], res_RTdiff[1], res_SA_RTdiff[1], res_SA_RTdiff_absRT[1]]
lst3 = [res_only_filtering[2], res_only_perc[2], res_SA[2], res_RTdiff[2], res_SA_RTdiff[2], res_SA_RTdiff_absRT[2]]
lst4 = [res_only_filtering[3], res_only_perc[3], res_SA[3], res_RTdiff[3], res_SA_RTdiff[3], res_SA_RTdiff_absRT[3]]
lst5 = [res_only_filtering[4], res_only_perc[4], res_SA[4], res_RTdiff[4], res_SA_RTdiff[4], res_SA_RTdiff_absRT[4]]
lst6 = [res_only_filtering[5], res_only_perc[5], res_SA[5], res_RTdiff[5], res_SA_RTdiff[5], res_SA_RTdiff_absRT[5]]
lst7 = [res_only_filtering[6], res_only_perc[6], res_SA[6], res_RTdiff[6], res_SA_RTdiff[6], res_SA_RTdiff_absRT[6]]
lst8 = [res_only_filtering[7], res_only_perc[7], res_SA[7], res_RTdiff[7], res_SA_RTdiff[7], res_SA_RTdiff_absRT[7]]
lst9 = [res_only_filtering[8], res_only_perc[8], res_SA[8], res_RTdiff[8], res_SA_RTdiff[8], res_SA_RTdiff_absRT[8]]
lst10 = [res_only_filtering[9], res_only_perc[9], res_SA[9], res_RTdiff[9], res_SA_RTdiff[9], res_SA_RTdiff_absRT[9]]
lst11 = [res_only_filtering[10], res_only_perc[10], res_SA[10], res_RTdiff[10], res_SA_RTdiff[10], res_SA_RTdiff_absRT[10]]
lst12 = [res_only_filtering[11], res_only_perc[11], res_SA[11], res_RTdiff[11], res_SA_RTdiff[11], res_SA_RTdiff_absRT[11]]


res_only_filtering_msgf = [21636,22661,20352,20179,21748,20460,20415,21692,21175,19514,21165,20122]
res_only_perc_msgf = [21641,22657,20362,20214,21736,20462,20412,21678,21200,19500,21171,20126]
res_SA_msgf = [21641,22657,20362,20214,21736,20462,20412,21678,21200,19500,21171,20126]
res_RTdiff_msgf = [21641,22770,20758,20340,22145,20866,20412,21678,21567,19500,21171,20126]
res_SA_RTdiff_msgf = [22096,23125,21113,20593,22475,21104,20859,22028,21771,20006,21463,20628]
res_SA_RTdiff_absRT_msgf = [22078,23067,21100,20585,22494,21103,20867,22021,21765,20033,21404,20601]
res_RTdiff_absRT_msgf = [21641,22657,20752,20214,22166,20886,20412,21678,21528,19500,21171,20126]


res_only_filtering_xtandem = [20551,21465,18785,18767,20605,19069,19217,20521,19746,18106,19909,18278]
res_only_perc_xtandem = [21123,22136,19491,19509,21057,19717,19807,21129,20473,18599,20610,18778]
res_SA_xtandem = [21662,22529,19888,19981,21671,20243,20187,21773,20958,19264,21191,19493]
res_RTdiff_xtandem = [22711,23776,21248,21163,22741,21270,21391,22656,22076,20375,22054,20913]
res_SA_RTdiff_xtandem = [23079,23952,21553,21443,23123,21606,21857,23111,22573,20717,22460,21235]
res_SA_RTdiff_absRT_xtandem = [23173,24055,21649,21451,23113,21790,21840,23192,22734,20733,22436,21366]
res_RTdiff_absRT_xtandem = [22868,23886,21368,21131,22765,21418,21459,22742,22252,20417,22067,21093]





print(lst1)

res_only_filtering_mean = sum(res_only_filtering) / len(res_only_filtering)
res_only_perc_mean = sum(res_only_perc) / len(res_only_perc)
res_SA_mean = sum(res_SA) / len(res_SA)
res_RTdiff_mean = sum(res_RTdiff) / len(res_RTdiff)
res_SA_RTdiff_mean = sum(res_SA_RTdiff) / len(res_SA_RTdiff)
res_SA_RTdiff_absRT_mean = sum(res_SA_RTdiff_absRT) / len(res_SA_RTdiff_absRT)


res_only_filtering_mean_top10 = sum(res_only_filtering_top10) / len(res_only_filtering_top10)
res_only_perc_mean_top10 = sum(res_only_perc_top10) / len(res_only_perc_top10)
res_SA_mean_top10 = sum(res_SA_top10) / len(res_SA_top10)
res_RTdiff_mean_top10 = sum(res_RTdiff_top10) / len(res_RTdiff_top10)
res_SA_RTdiff_mean_top10 = sum(res_SA_RTdiff_top10) / len(res_SA_RTdiff_top10)
res_SA_RTdiff_absRT_mean_top10 = sum(res_SA_RTdiff_absRT_top10) / len(res_SA_RTdiff_absRT_top10)


res_only_filtering_mean_msgf = sum(res_only_filtering_msgf) / len(res_only_filtering_msgf)
res_only_perc_mean_msgf = sum(res_only_perc_msgf) / len(res_only_perc_msgf)
res_SA_mean_msgf = sum(res_SA_msgf) / len(res_SA_msgf)
res_RTdiff_mean_msgf = sum(res_RTdiff_msgf) / len(res_RTdiff_msgf)
res_SA_RTdiff_mean_msgf = sum(res_SA_RTdiff_msgf) / len(res_SA_RTdiff_msgf)
res_SA_RTdiff_absRT_mean_msgf = sum(res_SA_RTdiff_absRT_msgf) / len(res_SA_RTdiff_absRT_msgf)


res_only_filtering_mean_xtandem = sum(res_only_filtering_xtandem) / len(res_only_filtering_xtandem)
res_only_perc_mean_xtandem = sum(res_only_perc_xtandem) / len(res_only_perc_xtandem)
res_SA_mean_xtandem = sum(res_SA_xtandem) / len(res_SA_xtandem)
res_RTdiff_mean_xtandem = sum(res_RTdiff_xtandem) / len(res_RTdiff_xtandem)
res_SA_RTdiff_mean_xtandem = sum(res_SA_RTdiff_xtandem) / len(res_SA_RTdiff_xtandem)
res_SA_RTdiff_absRT_mean_xtandem = sum(res_SA_RTdiff_absRT_xtandem) / len(res_SA_RTdiff_absRT_xtandem)


sse_mean_lst = [res_only_filtering_mean, res_only_perc_mean, res_SA_mean, res_RTdiff_mean, res_SA_RTdiff_mean, res_SA_RTdiff_absRT_mean]
msgf_mean_lst = [res_only_filtering_mean_msgf, res_only_perc_mean_msgf, res_SA_mean_msgf, res_RTdiff_mean_msgf, res_SA_RTdiff_mean_msgf, res_SA_RTdiff_absRT_mean_msgf]
xtandem_mean_lst = [res_only_filtering_mean_xtandem, res_only_perc_mean_xtandem, res_SA_mean_xtandem, res_RTdiff_mean_xtandem, res_SA_RTdiff_mean_xtandem, res_SA_RTdiff_absRT_mean_xtandem]


res_SA_sum = sum(res_SA)
res_RTdiff_sum = sum(res_RTdiff)
res_SA_RTdiff_sum = sum(res_SA_RTdiff)
res_SA_RTdiff_absRT_sum = sum(res_SA_RTdiff_absRT)



#hits_after_filtering_no_perc = [18750, 18800, 18809, 18810, 18853, 18874, 18876, 18856, 18810, 18774, 18694, 18535, 18481, 18349, 18232, 18140, 18071, 17963, 17773, 17601, 17186, 17096, 16938, 16707, 16547, 16195, 15833, 14964, 14852, 14733, 14589, 14493, 14424, 14385, 14249, 14211, 13750, 13789, 13827, 13656, 13201, 13119, 12816, 12551, 12180, 11870, 11209, 9244, 6453, 2619]
#hits_after_filtering_noiso = [19228,19263,19274,19280,19325,19263,19248,19327,19321,19292,19242,19211,19138,19117,19035,18988,18902,18735,18635,18542,18333,18221,18011,17871,17795,17617,17429,17288,17141,16934,16604,16547,16364,16157,16040,15904,15777,15690,15620,15656,15572,15506,15510,15336,15104,14625,14015,13562,12713,11540]
#hits_after_filtering = [19777,19846,19840,19900,19933,19941,19896,19882,19913,19837,19814,19692,19684,19601,19591,19503,19434,19378,19164,19121,18942,18807,18639,18487,18342,18285,18066,17734,17618,17371,17224,17019,16795,16641,16616,16510,16234,15937,15949,15943,15873,15830,15743,15592,15335,14802,14306,13643,12750,11727]

"""
x = filenames

plt.plot(x, res_only_filtering, label='Baseline')
#plt.plot(x, res_only_perc, label='Only Percolator')
#plt.plot(x, res_SA, label='SA')
#plt.plot(x, res_RTdiff, label='RT Difference')
#plt.plot(x, res_SA_RTdiff, label='SA + RT Difference')
plt.plot(x, res_SA_RTdiff_absRT, label='SA + RT Difference + predicted RT')
plt.legend()

#plt.xticks(range(0, 500, 20))


plt.xlabel('iPRG sample file')
plt.ylabel('# Hits after 1% FDR filtering')
plt.title('# Hits for the iPRG sample files')
#plt.xscale("log")
#plt.xlim(10, 500)
plt.ylim(16500, 21000)
plt.grid(True)
#plt.show()
plt.savefig('/vol/VMvol/data/PyCharm/basescript/deepLC/hits_deepLC_2curves_baseline.png')
"""


plt.xlabel('Integrated features\n')
plt.ylabel('Average # hits after 1% FDR filtering')
plt.title('Average # hits for the iPRG sample files\n(SSE, Top 10, 10 ppm prec. mass tol.)')
feats = ["Baseline", "Percolator", "+SA", "+$\delta$RT", "+SA+$\delta$RT", "+SA+$\delta$RT+|RT|"]
means = (res_only_filtering_mean_top10, res_only_perc_mean_top10, res_SA_mean_top10, res_RTdiff_mean_top10, res_SA_RTdiff_mean_top10, res_SA_RTdiff_absRT_mean_top10)
index = np.arange(len(feats))
#plt.legend(["Baseline", "Only Perc.", " + SA"], ['label1', 'label2', 'label3'])
plt.bar(index, means, 0.7)
plt.xticks(index, feats, fontsize=8)
plt.yticks(fontsize=9)
plt.ylim(17000, 22500)
plt.grid(True, axis='y')
#plt.tight_layout()
#plt.show()
plt.savefig('/vol/VMvol/data/PyCharm/basescript/deepLC/hits_avg_num_SSE_top10_sr.png')


"""
x = np.arange(5)
y1 = [34, 56, 12, 89, 67]
y2 = [12, 56, 78, 45, 90]
y3 = [14, 23, 45, 25, 89]
width = 0.2

# plot data in grouped manner of bar type
plt.bar(x - 0.2, y1, width, color='cyan')
plt.bar(x, y2, width, color='orange')
plt.bar(x + 0.2, y3, width, color='green')
plt.xticks(x, ['Team A', 'Team B', 'Team C', 'Team D', 'Team E'])
plt.xlabel("Teams")
plt.ylabel("Scores")
plt.legend(["Round 1", "Round 2", "Round 3"])
plt.show()
"""

"""
x = np.arange(6)
width = 0.2
feats = ["Baseline", "Percolator", "+SA", "+$\delta$RT", "+SA+$\delta$RT", "+SA+$\delta$RT+|RT|"]

plt.xlabel('Integrated features')
plt.ylabel('Average # hits after 1% FDR filtering')
plt.title('Average # hits for the iPRG sample files\n')

plt.bar(x - 0.2, sse_mean_lst, width)
plt.bar(x, msgf_mean_lst, width)
plt.bar(x + 0.2, xtandem_mean_lst, width)
index = np.arange(3)
plt.xticks(x, feats, fontsize=8)
plt.yticks(fontsize=9)
plt.ylim(17000, 22500)
plt.legend(["SSE", "MSGF+", "XTandem!"])
plt.grid(True, axis='y')

#plt.show()
plt.savefig('/vol/VMvol/data/PyCharm/basescript/deepLC/avg_hits_search_engines_grouped_sr.png')
"""

"""
x = np.arange(12)
width = 0.15

plt.xlabel('iPRG sample file')
plt.ylabel('# Hits after 1% FDR filtering')
plt.title('# Hits for the iPRG sample files\n(SSE, Top 10, 10 ppm prec. mass tol.)')
#feats = ["Baseline", "Percolator", " + SA", "+ RT Diff.", "+ Both", "+ abs. RT"]
feats = filenames
#means = (res_only_filtering_mean, res_only_perc_mean, res_SA_mean, res_RTdiff_mean, res_SA_RTdiff_mean, res_SA_RTdiff_absRT_mean)
index = np.arange(len(feats))
#plt.legend(["Baseline", "Only Perc.", " + SA"], ['label1', 'label2', 'label3'])
plt.bar(x - (width/2 + 2*width), res_only_filtering_top10, width)
plt.bar(x - (width/2 + width), res_only_perc_top10, width)
plt.bar(x - width/2, res_SA_top10, width)
plt.bar(x + width/2, res_RTdiff_top10, width)
plt.bar(x + (width/2 + width), res_SA_RTdiff_top10, width)
plt.bar(x + (width/2 + 2*width), res_SA_RTdiff_absRT_top10, width)
plt.xticks(index, feats)
plt.ylim(15000, 24250)
plt.legend(["Baseline", "Percolator", "+ SA", "+ $\delta$RT", "+ SA + $\delta$RT", "+ SA + $\delta$RT + |RT|"])
plt.grid(True, axis='y')
#plt.show()
plt.savefig('/vol/VMvol/data/PyCharm/basescript/deepLC/hits_grouped_SSE_top10_sr.png')
"""

"""
x = np.arange(6)
width = 0.75

plt.xlabel('iPRG sample file')
plt.ylabel('# Hits after 1% FDR filtering')
plt.title('# Hits for the iPRG sample files')
feats = ["Baseline", "Percolator", " + SA", "+ RT Diff.", "+ Both", "+ abs. RT"]
#feats = filenames
means = (res_only_filtering_mean, res_only_perc_mean, res_SA_mean, res_RTdiff_mean, res_SA_RTdiff_mean, res_SA_RTdiff_absRT_mean)
index = np.arange(len(feats))
#plt.legend(["Baseline", "Only Perc.", " + SA"], ['label1', 'label2', 'label3'])
plt.bar(x - (width/2 + 5*width), lst1, width)
plt.bar(x - (width/2 + 4*width), lst2, width)
plt.bar(x - (width/2 + 3*width), lst3, width)
plt.bar(x - (width/2 + 2*width), lst4, width)
plt.bar(x - (width/2 + width), lst5, width)
plt.bar(x - width/2, lst6, width)
plt.bar(x + width/2, lst7, width)
plt.bar(x + (width/2 + width), lst8, width)
plt.bar(x + (width/2 + 2*width), lst9, width)
plt.bar(x + (width/2 + 3*width), lst10, width)
plt.bar(x + (width/2 + 4*width), lst11, width)
plt.bar(x + (width/2 + 5*width), lst12, width)
plt.xticks(index, feats)
plt.ylim(16500, 21250)
#plt.legend(["Baseline", "Percolator", "+ SA", "+ RT Diff.", "+ Both", "+ abs. RT"])
plt.grid(True, axis='y')
plt.show()
#plt.savefig('/vol/VMvol/data/PyCharm/basescript/deepLC/hits_deepLC_grouped.png')
"""


"""
x = np.arange(12)
width = 0.2

plt.xlabel('iPRG sample file')
plt.ylabel('# Hits after 1% FDR filtering')
plt.title('# Hits for the iPRG sample files\n(SSE, Top 10, 10 ppm prec. mass tol.)')
#feats = ["Baseline", "Percolator", " + SA", "+ RT Diff.", "+ Both", "+ abs. RT"]
feats = filenames
#means = (res_only_filtering_mean, res_only_perc_mean, res_SA_mean, res_RTdiff_mean, res_SA_RTdiff_mean, res_SA_RTdiff_absRT_mean)
index = np.arange(len(feats))
#plt.legend(["Baseline", "Only Perc.", " + SA"], ['label1', 'label2', 'label3'])
plt.bar(x - (width/2 + width), res_only_filtering_top10, width)
plt.bar(x - width/2, res_only_perc_top10, width)
#plt.bar(x - width/2, res_SA, width)
plt.bar(x + width/2, res_RTdiff_top10, width)
#plt.bar(x + (width/2 + width), res_SA_RTdiff, width)
plt.bar(x + (width/2 + width), res_RTdiff_absRT_top10, width)
plt.xticks(index, feats)
plt.ylim(15000, 24250)
plt.legend(["Baseline", "Percolator", "+ $\delta$RT", "+ $\delta$RT + |RT|"])
plt.grid(True, axis='y')
#plt.show()
plt.savefig('/vol/VMvol/data/PyCharm/basescript/deepLC/hits_grouped_onlyRT_SSE_top10_sr.png')
"""
