###This file includes codes for performance evaluation###

###Calculate ARI and V-measure in python###
from sklearn import metrics

def ARI(pre_label1, pre_label2):
    return(metrics.adjusted_rand_score(pre_label1, pre_label2))

def V_measure(pre_label1, pre_label2):
    return(metrics.v_measure_score(pre_label1, pre_label2))
