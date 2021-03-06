delfile <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/Elevation/Fst_less2500_vsmore2500.txt",header = T, sep="\t")

head(delfile)
quantile(delfile$FST,probs = seq(0,1,0.005),na.rm = T)

SUBSET <- delfile[(delfile[4] >= 0.04180882),]
nrow (SUBSET)

write.table(SUBSET, file = "/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/Elevation/Outliers/Fst_less2500_vsmore2500_outliers.txt",row.names=F,col.names=F,sep="\t",quote=F)

#get the outliers for the file of Fst_less3000_vsmore3000.txt:

elevation_3000 <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/Elevation/Fst_less3000_vsmore3000_no_NA.txt",header = T, sep="\t")

head(elevation_3000)
quantile(elevation_3000$FST,probs = seq(0,1,0.005),na.rm = T)

# 0.0%          0.5%          1.0%          1.5%          2.0%          2.5%          3.0%          3.5%          4.0%          4.5%          5.0%          5.5%          6.0% 
#-0.0130300590 -0.0100346060 -0.0099755465 -0.0099010243 -0.0097949811 -0.0096503637 -0.0095232039 -0.0093099973 -0.0092392882 -0.0090612959 -0.0088409685 -0.0085906271 -0.0084835267 
#6.5%          7.0%          7.5%          8.0%          8.5%          9.0%          9.5%         10.0%         10.5%         11.0%         11.5%         12.0%         12.5% 
#-0.0081687055 -0.0079540535 -0.0077939893 -0.0075132247 -0.0071417978 -0.0069246553 -0.0066556224 -0.0063716840 -0.0061595415 -0.0058776655 -0.0055358580 -0.0053012623 -0.0048828407 
#13.0%         13.5%         14.0%         14.5%         15.0%         15.5%         16.0%         16.5%         17.0%         17.5%         18.0%         18.5%         19.0% 
#-0.0046388000 -0.0041012326 -0.0038309449 -0.0033943659 -0.0032348099 -0.0026500518 -0.0021723546 -0.0017730163 -0.0014585425 -0.0010063770 -0.0006990090 -0.0002329290  0.0000977159 
#19.5%         20.0%         20.5%         21.0%         21.5%         22.0%         22.5%         23.0%         23.5%         24.0%         24.5%         25.0%         25.5% 
#0.0005284770  0.0009209356  0.0013047319  0.0018611326  0.0022262127  0.0027218589  0.0030738201  0.0035837850  0.0039720712  0.0044954890  0.0049796495  0.0052532370  0.0057781099 
#26.0%         26.5%         27.0%         27.5%         28.0%         28.5%         29.0%         29.5%         30.0%         30.5%         31.0%         31.5%         32.0% 
#0.0062337306  0.0066839421  0.0072958272  0.0078909243  0.0084425860  0.0090108270  0.0097158004  0.0101591538  0.0108221890  0.0118530954  0.0124904067  0.0128886284  0.0135886076 
#32.5%         33.0%         33.5%         34.0%         34.5%         35.0%         35.5%         36.0%         36.5%         37.0%         37.5%         38.0%         38.5% 
#0.0141863732  0.0149544688  0.0155542680  0.0161591663  0.0172099760  0.0177141869  0.0184452493  0.0191530836  0.0198561311  0.0205369700  0.0214521250  0.0221495854  0.0229856802 
#39.0%         39.5%         40.0%         40.5%         41.0%         41.5%         42.0%         42.5%         43.0%         43.5%         44.0%         44.5%         45.0% 
#0.0239771939  0.0247085899  0.0257195274  0.0264503415  0.0270037160  0.0276944049  0.0284844566  0.0293159751  0.0300946365  0.0310288141  0.0319379194  0.0329992831  0.0335959222 
#45.5%         46.0%         46.5%         47.0%         47.5%         48.0%         48.5%         49.0%         49.5%         50.0%         50.5%         51.0%         51.5% 
#0.0343287868  0.0353004913  0.0360400959  0.0369384162  0.0376819607  0.0386147086  0.0393405871  0.0401366560  0.0412269720  0.0420629480  0.0428246415  0.0435464859  0.0440713741 
#52.0%         52.5%         53.0%         53.5%         54.0%         54.5%         55.0%         55.5%         56.0%         56.5%         57.0%         57.5%         58.0% 
#0.0450922462  0.0459604718  0.0465158734  0.0471163634  0.0479582850  0.0486005412  0.0494951383  0.0506435498  0.0513632400  0.0528204073  0.0534347192  0.0542662548  0.0555150830 
#58.5%         59.0%         59.5%         60.0%         60.5%         61.0%         61.5%         62.0%         62.5%         63.0%         63.5%         64.0%         64.5% 
#0.0564385257  0.0572770762  0.0581707557  0.0595237158  0.0605467488  0.0614704126  0.0627142807  0.0637056995  0.0646046810  0.0656959882  0.0672761743  0.0689622949  0.0701269230 
#65.0%         65.5%         66.0%         66.5%         67.0%         67.5%         68.0%         68.5%         69.0%         69.5%         70.0%         70.5%         71.0% 
#0.0712255323  0.0723817106  0.0733201219  0.0740795241  0.0756540780  0.0770057347  0.0786325203  0.0801886867  0.0814408573  0.0827478103  0.0841749738  0.0856662106  0.0871116353 
#71.5%         72.0%         72.5%         73.0%         73.5%         74.0%         74.5%         75.0%         75.5%         76.0%         76.5%         77.0%         77.5% 
#0.0882751512  0.0897654456  0.0912770196  0.0926138335  0.0945085375  0.0964641950  0.0982006559  0.1003334060  0.1022608334  0.1038552831  0.1060215460  0.1086960073  0.1108619960 
#78.0%         78.5%         79.0%         79.5%         80.0%         80.5%         81.0%         81.5%         82.0%         82.5%         83.0%         83.5%         84.0% 
#0.1137295893  0.1154951588  0.1174984995  0.1195801309  0.1217891972  0.1242849741  0.1272730099  0.1304298752  0.1326311092  0.1347622658  0.1379083063  0.1403514926  0.1443256012 
#84.5%         85.0%         85.5%         86.0%         86.5%         87.0%         87.5%         88.0%         88.5%         89.0%         89.5%         90.0%         90.5% 
#0.1488119224  0.1518883241  0.1546586409  0.1579605094  0.1611069763  0.1636060047  0.1681814233  0.1723862926  0.1756082640  0.1788384189  0.1823164973  0.1855908612  0.1912696631 
#91.0%         91.5%         92.0%         92.5%         93.0%         93.5%         94.0%         94.5%         95.0%         95.5%         96.0%         96.5%         97.0% 
#0.1974153155  0.2022329140  0.2088395602  0.2163907403  0.2217591437  0.2275819771  0.2361519163  0.2428144328  0.2498974176  0.2586623386  0.2653058776  0.2765785923  0.2868774253 
#97.5%         98.0%         98.5%         99.0%         99.5%        100.0% 
#  0.3014466887  0.3259427894  0.3501763844  0.3810938052  0.4177053646  0.6492715560 
SUBSET <- elevation_3000[(elevation_3000[4] >= 0.3810938052),]
nrow (SUBSET)

write.table(SUBSET, file = "/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/Elevation/Outliers/Fst_less3000_vsmore3000_outliers_99th.txt",row.names=F,col.names=F,sep="\t",quote=F)


#get the outliers for the file of Fst_wild_range30_40_vs_higherLat40.txt:
wild_higher40 <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/Latitude/Fst_wild_range30_40_vs_higherLat40_no_NA.txt",header = T, sep="\t")

head(wild_higher40)
quantile(wild_higher40$FST,probs = seq(0,1,0.005),na.rm = T)

#  0.0%          0.5%          1.0%          1.5%          2.0%          2.5%          3.0%          3.5%          4.0%          4.5%          5.0%          5.5%          6.0% 
#-0.0036138000 -0.0033090500 -0.0032921518 -0.0032696414 -0.0032385683 -0.0032089633 -0.0031768800 -0.0031432378 -0.0030840812 -0.0030196073 -0.0029419718 -0.0028553458 -0.0027473597 
#6.5%          7.0%          7.5%          8.0%          8.5%          9.0%          9.5%         10.0%         10.5%         11.0%         11.5%         12.0%         12.5% 
#-0.0026495073 -0.0025346154 -0.0024175345 -0.0023102063 -0.0021984502 -0.0021177693 -0.0019914290 -0.0018950133 -0.0017842035 -0.0016177104 -0.0014958628 -0.0013509570 -0.0011572414 
#13.0%         13.5%         14.0%         14.5%         15.0%         15.5%         16.0%         16.5%         17.0%         17.5%         18.0%         18.5%         19.0% 
#-0.0009697043 -0.0008377296 -0.0006589150 -0.0005665999 -0.0003937901 -0.0002068069 -0.0000303280  0.0002326203  0.0004265556  0.0006861620  0.0008507559  0.0011930704  0.0015750422 
#19.5%         20.0%         20.5%         21.0%         21.5%         22.0%         22.5%         23.0%         23.5%         24.0%         24.5%         25.0%         25.5% 
#0.0017876365  0.0020273200  0.0022450791  0.0025266861  0.0027265930  0.0030420928  0.0033603484  0.0036480283  0.0040119685  0.0043238484  0.0046436148  0.0049426160  0.0052394468 
#26.0%         26.5%         27.0%         27.5%         28.0%         28.5%         29.0%         29.5%         30.0%         30.5%         31.0%         31.5%         32.0% 
#0.0056384335  0.0061903959  0.0066381098  0.0069295433  0.0072717572  0.0076724412  0.0080501720  0.0083793010  0.0088651630  0.0092060685  0.0097394274  0.0100980610  0.0104851402 
#32.5%         33.0%         33.5%         34.0%         34.5%         35.0%         35.5%         36.0%         36.5%         37.0%         37.5%         38.0%         38.5% 
#0.0109610406  0.0114412967  0.0118918701  0.0124315825  0.0127667051  0.0131984624  0.0137460350  0.0141215208  0.0147627597  0.0153591961  0.0158238043  0.0165083354  0.0170719520 
#39.0%         39.5%         40.0%         40.5%         41.0%         41.5%         42.0%         42.5%         43.0%         43.5%         44.0%         44.5%         45.0% 
#0.0177521494  0.0184173940  0.0190512766  0.0195458440  0.0203093704  0.0209447417  0.0214810052  0.0222495509  0.0230828835  0.0236355828  0.0241775864  0.0247631379  0.0253219844 
#45.5%         46.0%         46.5%         47.0%         47.5%         48.0%         48.5%         49.0%         49.5%         50.0%         50.5%         51.0%         51.5% 
#0.0258049679  0.0268515050  0.0274536734  0.0282646110  0.0289302386  0.0295812307  0.0300079294  0.0307187962  0.0313956018  0.0319556585  0.0324657874  0.0329374167  0.0338581774 
#52.0%         52.5%         53.0%         53.5%         54.0%         54.5%         55.0%         55.5%         56.0%         56.5%         57.0%         57.5%         58.0% 
#0.0346448320  0.0351394815  0.0359267249  0.0367821270  0.0377979612  0.0385966139  0.0394611220  0.0401140283  0.0410386532  0.0418106656  0.0428296021  0.0438178905  0.0445883924 
#58.5%         59.0%         59.5%         60.0%         60.5%         61.0%         61.5%         62.0%         62.5%         63.0%         63.5%         64.0%         64.5% 
#0.0451815901  0.0462762683  0.0469725803  0.0480152326  0.0488708176  0.0497588149  0.0507926039  0.0519644298  0.0529671256  0.0540288974  0.0551798387  0.0560765310  0.0571110855 
#65.0%         65.5%         66.0%         66.5%         67.0%         67.5%         68.0%         68.5%         69.0%         69.5%         70.0%         70.5%         71.0% 
#0.0582637229  0.0593424423  0.0606592445  0.0616528978  0.0628722134  0.0644345427  0.0661209497  0.0675515169  0.0683973822  0.0699640849  0.0714351406  0.0728049605  0.0736980062 
#71.5%         72.0%         72.5%         73.0%         73.5%         74.0%         74.5%         75.0%         75.5%         76.0%         76.5%         77.0%         77.5% 
#0.0753157235  0.0761845109  0.0777833276  0.0793990122  0.0809346051  0.0823292738  0.0840634203  0.0856684680  0.0876606283  0.0890116724  0.0905225500  0.0919994374  0.0937157204 
#78.0%         78.5%         79.0%         79.5%         80.0%         80.5%         81.0%         81.5%         82.0%         82.5%         83.0%         83.5%         84.0% 
#0.0950372503  0.0968552597  0.0993637961  0.1014369596  0.1032984624  0.1055254086  0.1073788392  0.1096427282  0.1127403960  0.1146686536  0.1174121349  0.1202899123  0.1221943657 
#84.5%         85.0%         85.5%         86.0%         86.5%         87.0%         87.5%         88.0%         88.5%         89.0%         89.5%         90.0%         90.5% 
#0.1240694696  0.1278987249  0.1313019882  0.1338469970  0.1373639373  0.1393098857  0.1419366344  0.1444600450  0.1486595703  0.1526939671  0.1572727994  0.1616525565  0.1653529279 
#91.0%         91.5%         92.0%         92.5%         93.0%         93.5%         94.0%         94.5%         95.0%         95.5%         96.0%         96.5%         97.0% 
#0.1696226773  0.1759195386  0.1787973178  0.1820943338  0.1854240665  0.1923093162  0.1958941808  0.2017946429  0.2088322589  0.2153318809  0.2249295036  0.2361146704  0.2461237740 
#97.5%         98.0%         98.5%         99.0%         99.5%        100.0% 
#  0.2601220106  0.2763243143  0.3022058331  0.3298958971  0.3723095788  0.4839167440 


SUBSET <- wild_higher40[(wild_higher40[4] >= 0.3298958971 ),]
nrow (SUBSET)

write.table(SUBSET, file = "/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/Elevation/Outliers/Fst_wild_range30_40_vs_higherLat40_outliers_99th.txt",row.names=F,col.names=F,sep="\t",quote=F)


##get the outliers for the file of Fst_wild_range30_40_vs_lowerLat30.txt:
wild_lower30 <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/Latitude/Fst_moreoreq30to40vsless30_no_NA.txt",header = T, sep="\t")

head(wild_lower30)
quantile(wild_lower30$FST,probs = seq(0,1,0.005),na.rm = F)

#0.0%          0.5%          1.0%          1.5%          2.0%          2.5%          3.0%          3.5%          4.0%          4.5%          5.0%          5.5%          6.0% 
#-0.0063750240 -0.0044311020 -0.0044220765 -0.0043765992 -0.0043434263 -0.0042982298 -0.0042577416 -0.0041968102 -0.0041146770 -0.0040340587 -0.0039705430 -0.0038997360 -0.0037809421 
#6.5%          7.0%          7.5%          8.0%          8.5%          9.0%          9.5%         10.0%         10.5%         11.0%         11.5%         12.0%         12.5% 
#-0.0036201591 -0.0034636830 -0.0033191980 -0.0032550411 -0.0031119829 -0.0028874877 -0.0027300174 -0.0025128774 -0.0022479225 -0.0022172570 -0.0021770510 -0.0019134327 -0.0016982158 
#13.0%         13.5%         14.0%         14.5%         15.0%         15.5%         16.0%         16.5%         17.0%         17.5%         18.0%         18.5%         19.0% 
#-0.0013462652 -0.0010728600 -0.0007476894 -0.0004916937 -0.0002603840 -0.0000043000  0.0001802064  0.0005135499  0.0009185630  0.0012225560  0.0015311098  0.0019143832  0.0022210257 
#19.5%         20.0%         20.5%         21.0%         21.5%         22.0%         22.5%         23.0%         23.5%         24.0%         24.5%         25.0%         25.5% 
#0.0024018320  0.0025696196  0.0028540540  0.0032262718  0.0035770475  0.0040160680  0.0042409055  0.0044148306  0.0045013716  0.0048769522  0.0053288801  0.0057739070  0.0061730164 
#26.0%         26.5%         27.0%         27.5%         28.0%         28.5%         29.0%         29.5%         30.0%         30.5%         31.0%         31.5%         32.0% 
#0.0066296353  0.0070500251  0.0075144958  0.0077426260  0.0081678484  0.0086198010  0.0089541702  0.0094204270  0.0098888836  0.0104691485  0.0109583583  0.0113977143  0.0119062707 
#32.5%         33.0%         33.5%         34.0%         34.5%         35.0%         35.5%         36.0%         36.5%         37.0%         37.5%         38.0%         38.5% 
#0.0123715815  0.0129292020  0.0132992387  0.0137998384  0.0142799173  0.0148503766  0.0154617432  0.0158013778  0.0163491152  0.0170297210  0.0175513675  0.0181709218  0.0188991280 
#39.0%         39.5%         40.0%         40.5%         41.0%         41.5%         42.0%         42.5%         43.0%         43.5%         44.0%         44.5%         45.0% 
#0.0193525687  0.0199307129  0.0205843882  0.0211951141  0.0220091710  0.0225438177  0.0233511901  0.0240011400  0.0245846478  0.0251917576  0.0258353572  0.0266495988  0.0275130024 
#45.5%         46.0%         46.5%         47.0%         47.5%         48.0%         48.5%         49.0%         49.5%         50.0%         50.5%         51.0%         51.5% 
#0.0283085540  0.0290771639  0.0300221011  0.0306598575  0.0314548250  0.0321598020  0.0330165553  0.0338629226  0.0345159576  0.0352387230  0.0360767106  0.0370089147  0.0377523650 
#52.0%         52.5%         53.0%         53.5%         54.0%         54.5%         55.0%         55.5%         56.0%         56.5%         57.0%         57.5%         58.0% 
#0.0382402106  0.0391176669  0.0397526187  0.0404057725  0.0410302441  0.0421178328  0.0431366535  0.0442238168  0.0452164796  0.0461343580  0.0471735722  0.0484740620  0.0494144500 
#58.5%         59.0%         59.5%         60.0%         60.5%         61.0%         61.5%         62.0%         62.5%         63.0%         63.5%         64.0%         64.5% 
#0.0508262690  0.0521508910  0.0530136156  0.0542085684  0.0549661279  0.0563538382  0.0578908454  0.0588173733  0.0598396578  0.0610561827  0.0620260713  0.0631821512  0.0645792152 
#65.0%         65.5%         66.0%         66.5%         67.0%         67.5%         68.0%         68.5%         69.0%         69.5%         70.0%         70.5%         71.0% 
#0.0654353382  0.0669366454  0.0678581674  0.0693661389  0.0705675364  0.0718802302  0.0734009316  0.0744483604  0.0759151026  0.0770131145  0.0785610278  0.0801388852  0.0811861909 
#71.5%         72.0%         72.5%         73.0%         73.5%         74.0%         74.5%         75.0%         75.5%         76.0%         76.5%         77.0%         77.5% 
#0.0825144274  0.0841152836  0.0856204633  0.0873035761  0.0889430421  0.0898676487  0.0914178610  0.0932593960  0.0950764214  0.0966678328  0.0987767060  0.1010989490  0.1036793204 
#78.0%         78.5%         79.0%         79.5%         80.0%         80.5%         81.0%         81.5%         82.0%         82.5%         83.0%         83.5%         84.0% 
#0.1058611611  0.1080421993  0.1097860302  0.1116862472  0.1141332182  0.1174586745  0.1202227934  0.1230349672  0.1262053512  0.1292528096  0.1318092640  0.1338972641  0.1361343866 
#84.5%         85.0%         85.5%         86.0%         86.5%         87.0%         87.5%         88.0%         88.5%         89.0%         89.5%         90.0%         90.5% 
#0.1397144133  0.1429461582  0.1464194589  0.1507324157  0.1544008369  0.1590852499  0.1635552335  0.1674220672  0.1713721491  0.1762075842  0.1815004539  0.1867161780  0.1913561537 
#91.0%         91.5%         92.0%         92.5%         93.0%         93.5%         94.0%         94.5%         95.0%         95.5%         96.0%         96.5%         97.0% 
#0.1977108360  0.2027472016  0.2087146348  0.2130116513  0.2189672415  0.2263475382  0.2328214999  0.2392980908  0.2457957113  0.2568265407  0.2666874315  0.2797138686  0.2945335264 
#97.5%         98.0%         98.5%         99.0%         99.5%        100.0% 
#  0.3124423020  0.3371896812  0.3570533716  0.3874640773  0.4398500371  0.5171261540 

SUBSET <- wild_lower30[(wild_lower30[4] >= 0.3874640773 ),]
nrow (SUBSET)

write.table(SUBSET, file = "/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/Elevation/Outliers/Fst_moreoreq30to40vsless30_no_NA_outliers_99th.txt",row.names=F,col.names=F,sep="\t",quote=F)


##get the outliers for the file of Fst_wild_range30_40_vs_lowerLat30.txt:
GH <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/GrowthHabit/Fst_average_80samp_10iter_no_NA_physPos.txt",header = T, sep="\t")

head(GH)
quantile(GH$FST,probs = seq(0,1,0.005),na.rm = F)

#0.0%          0.5%          1.0%          1.5%          2.0%          2.5%          3.0%          3.5%          4.0%          4.5%          5.0%          5.5%          6.0% 
#-9.182963e-03 -7.730634e-03 -7.226883e-03 -6.761754e-03 -6.439952e-03 -6.169286e-03 -5.979320e-03 -5.760352e-03 -5.533124e-03 -5.312634e-03 -5.170848e-03 -4.990239e-03 -4.849380e-03 
#6.5%          7.0%          7.5%          8.0%          8.5%          9.0%          9.5%         10.0%         10.5%         11.0%         11.5%         12.0%         12.5% 
#-4.710378e-03 -4.595951e-03 -4.431939e-03 -4.275824e-03 -4.069258e-03 -3.918106e-03 -3.711077e-03 -3.526581e-03 -3.353930e-03 -3.235149e-03 -3.093287e-03 -2.909445e-03 -2.711051e-03 
#13.0%         13.5%         14.0%         14.5%         15.0%         15.5%         16.0%         16.5%         17.0%         17.5%         18.0%         18.5%         19.0% 
#-2.533163e-03 -2.461847e-03 -2.330017e-03 -2.185068e-03 -2.055024e-03 -1.847266e-03 -1.640193e-03 -1.467777e-03 -1.279962e-03 -1.095410e-03 -9.054014e-04 -7.281770e-04 -4.649428e-04 
#19.5%         20.0%         20.5%         21.0%         21.5%         22.0%         22.5%         23.0%         23.5%         24.0%         24.5%         25.0%         25.5% 
#-3.046988e-04 -4.058513e-05  1.710838e-04  3.635607e-04  5.271212e-04  7.780613e-04  1.023456e-03  1.280130e-03  1.552647e-03  1.763626e-03  2.046700e-03  2.270904e-03  2.542854e-03 
#26.0%         26.5%         27.0%         27.5%         28.0%         28.5%         29.0%         29.5%         30.0%         30.5%         31.0%         31.5%         32.0% 
#2.837213e-03  3.071307e-03  3.326594e-03  3.618728e-03  3.928155e-03  4.203384e-03  4.412132e-03  4.689928e-03  4.883016e-03  5.198140e-03  5.426610e-03  5.657052e-03  5.842833e-03 
#32.5%         33.0%         33.5%         34.0%         34.5%         35.0%         35.5%         36.0%         36.5%         37.0%         37.5%         38.0%         38.5% 
#6.089338e-03  6.442204e-03  6.715561e-03  7.065911e-03  7.467785e-03  7.891869e-03  8.346929e-03  8.691760e-03  9.148819e-03  9.513031e-03  9.961435e-03  1.048597e-02  1.089293e-02 
#39.0%         39.5%         40.0%         40.5%         41.0%         41.5%         42.0%         42.5%         43.0%         43.5%         44.0%         44.5%         45.0% 
#1.125310e-02  1.164922e-02  1.202275e-02  1.247055e-02  1.294123e-02  1.328100e-02  1.364336e-02  1.415981e-02  1.469342e-02  1.537238e-02  1.594946e-02  1.641704e-02  1.688580e-02 
#45.5%         46.0%         46.5%         47.0%         47.5%         48.0%         48.5%         49.0%         49.5%         50.0%         50.5%         51.0%         51.5% 
#1.751926e-02  1.813856e-02  1.878909e-02  1.939420e-02  1.993015e-02  2.038162e-02  2.094329e-02  2.148256e-02  2.216059e-02  2.258770e-02  2.311996e-02  2.374989e-02  2.449938e-02 
#52.0%         52.5%         53.0%         53.5%         54.0%         54.5%         55.0%         55.5%         56.0%         56.5%         57.0%         57.5%         58.0% 
#2.511098e-02  2.577927e-02  2.623960e-02  2.680764e-02  2.735521e-02  2.796117e-02  2.870227e-02  2.928226e-02  2.992376e-02  3.066001e-02  3.136646e-02  3.221110e-02  3.277277e-02 
#58.5%         59.0%         59.5%         60.0%         60.5%         61.0%         61.5%         62.0%         62.5%         63.0%         63.5%         64.0%         64.5% 
#3.326399e-02  3.403015e-02  3.492031e-02  3.570809e-02  3.639839e-02  3.721852e-02  3.784385e-02  3.853972e-02  3.970368e-02  4.042782e-02  4.148830e-02  4.211170e-02  4.324635e-02 
#65.0%         65.5%         66.0%         66.5%         67.0%         67.5%         68.0%         68.5%         69.0%         69.5%         70.0%         70.5%         71.0% 
#4.426725e-02  4.550631e-02  4.646580e-02  4.710425e-02  4.832735e-02  4.888370e-02  4.999941e-02  5.096335e-02  5.173371e-02  5.246681e-02  5.351052e-02  5.434644e-02  5.553314e-02 
#71.5%         72.0%         72.5%         73.0%         73.5%         74.0%         74.5%         75.0%         75.5%         76.0%         76.5%         77.0%         77.5% 
#5.657368e-02  5.784089e-02  5.899013e-02  5.985312e-02  6.062005e-02  6.169101e-02  6.296962e-02  6.422293e-02  6.604205e-02  6.752395e-02  6.895676e-02  7.020749e-02  7.155782e-02 
#78.0%         78.5%         79.0%         79.5%         80.0%         80.5%         81.0%         81.5%         82.0%         82.5%         83.0%         83.5%         84.0% 
#7.257619e-02  7.400993e-02  7.530337e-02  7.669865e-02  7.802632e-02  7.977014e-02  8.137679e-02  8.292953e-02  8.496618e-02  8.722806e-02  8.963085e-02  9.217496e-02  9.366112e-02 
#84.5%         85.0%         85.5%         86.0%         86.5%         87.0%         87.5%         88.0%         88.5%         89.0%         89.5%         90.0%         90.5% 
#9.490396e-02  9.694021e-02  9.866482e-02  1.013684e-01  1.031796e-01  1.051328e-01  1.072812e-01  1.090584e-01  1.118687e-01  1.137520e-01  1.166277e-01  1.193818e-01  1.231024e-01 
#91.0%         91.5%         92.0%         92.5%         93.0%         93.5%         94.0%         94.5%         95.0%         95.5%         96.0%         96.5%         97.0% 
#1.264955e-01  1.305371e-01  1.345630e-01  1.384751e-01  1.435495e-01  1.503444e-01  1.548805e-01  1.599707e-01  1.638698e-01  1.687138e-01  1.752450e-01  1.836968e-01  1.893786e-01 
#97.5%         98.0%         98.5%         99.0%         99.5%        100.0% 
#  1.991411e-01  2.087320e-01  2.259717e-01  2.542977e-01  2.898167e-01  3.948482e-01 

SUBSET <- GH[(GH[4] >= 0.2542977 ),]
nrow (SUBSET)

write.table(SUBSET, file = "/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/GrowthHabit/Outliers/outlier_Fst_SPRING_vs_WINTER_80samples_average_noNA_99th.txt",row.names=F,col.names=F,sep="\t",quote=F)
