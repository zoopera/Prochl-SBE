#!/usr/bin/env python

import os
import sys

# input a file and return the index (starting from 0) of a list of titles
def get_index(ffile, requested_title):
    # load the first line (title line)
    titles = []
    with open(ffile) as f:
	titles = f.readline()[:-1].split("\t")
    # find the index for each requested title
    ret_index = []
    for title in requested_title:
	i = 0
	while i<len(titles) and titles[i]!=title:
	    i += 1
	if i<len(titles):
	    ret_index.append(i)
	else:
	    ret_index.append(10000)	# if title not found, a faulty value of 10000 will be assigned
    return ret_index


if __name__ == "__main__":

  input = sys.argv[1]
  group = sys.argv[2]  # "Alphaproteobacteria"

  # get index of titles
  [ind_class, ind_clade, ind_is_target, ind_cmp_clade] = get_index(input, ["classification", "clade",  "is_target", "cmp_clade"])
  [ind_HON3_avg, ind_HON3_se, ind_HON3_signtest, ind_HON3_ttest] = get_index(input, ["HON3_avg", "HON3_se", "HON3_sign_test", "HON3_ttest"])
  [ind_HON2_avg, ind_HON2_se, ind_HON2_signtest, ind_HON2_ttest] = get_index(input, ["HON2_avg", "HON2_se", "HON2_sign_test", "HON2_ttest"])
  [ind_HON0_avg, ind_HON0_se, ind_HON0_signtest, ind_HON0_ttest] = get_index(input, ["HON0_avg", "HON0_se", "HON0_sign_test", "HON0_ttest"])

  buff = "\t\t\tRCCulculator vs3 (GC-corrected based on codon frequency)\t\t\t\t"
  buff += "RCCulculator vs2 (GC-corrected based on AA frequency)\t\t\t\t"
  buff += "RCCulculator vs1 (uncorrected)\t\t\t\n"
  buff += "Group\tTarget Clade\tControl Clade\tdRdC (target)\tdRdC (control)\tp-value 1\tp-value 2\t"
  buff += "dRdC (target)\tdRdC (control)\tp-value 1\tp-value 2\t"
  buff += "dRdC (target)\tdRdC (control)\tp-value 1\tp-value 2\n"

  for classi in ["charge", "MY"]:

    # load dRdC values of target
    dict_target_dRdC = {} # indexed by control clade
    with open(input) as f:
        line = f.readline()
	while line:
	    arr = line[:-1].split("\t")
	    if arr[ind_class]==classi and arr[ind_is_target]=="1" and  arr[ind_cmp_clade]!="NA" :
		dict_target_dRdC[arr[ind_cmp_clade]]={}
		dict_target_dRdC[arr[ind_cmp_clade]]["HON3"] = "%s&%s" % (arr[ind_HON3_avg], arr[ind_HON3_se])
		dict_target_dRdC[arr[ind_cmp_clade]]["HON2"] = "%s&%s" % (arr[ind_HON2_avg], arr[ind_HON2_se])
	  	dict_target_dRdC[arr[ind_cmp_clade]]["HON0"] = "%s&%s" % (arr[ind_HON0_avg], arr[ind_HON0_se])
	    line = f.readline()

    # load control clade values and write to file
    buff += "%s (%s)" % (group, classi)
    with open(input) as f:
	line = f.readline()
	while line:
	    arr = line[:-1].split("\t")
	    if arr[ind_class]==classi and arr[ind_is_target]=="0":
		target = arr[ind_cmp_clade]
		control = arr[ind_clade]
		buff += "\t%s\t%s\t" % (target, control)
		buff += "%s\t%s&%s\t%s\t%s\t" % (dict_target_dRdC[control]["HON3"], arr[ind_HON3_avg], arr[ind_HON3_se], arr[ind_HON3_signtest], arr[ind_HON3_ttest])
		buff += "%s\t%s&%s\t%s\t%s\t" % (dict_target_dRdC[control]["HON2"], arr[ind_HON2_avg], arr[ind_HON2_se], arr[ind_HON2_signtest], arr[ind_HON2_ttest])
		buff += "%s\t%s&%s\t%s\t%s\n" % (dict_target_dRdC[control]["HON0"], arr[ind_HON0_avg], arr[ind_HON0_se], arr[ind_HON0_signtest], arr[ind_HON0_ttest])
	    line = f.readline()


  with open("%s.formatted" % (input), "w") as f:
    f.write(buff)

