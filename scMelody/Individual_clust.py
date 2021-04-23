from sklearn.cluster import SpectralClustering
from sklearn import metrics

def sc_pre(sm, C):
    sc = SpectralClustering(n_clusters=C, n_neighbors=2*C, affinity='precomputed',n_init=100,assign_labels='kmeans')
    out = sc.fit_predict(sm)
    return out

def sil_score(sm, pre_lable):
    score1 = metrics.silhouette_score(1-sm, out, metric='precomputed')
    return np.array([score1])

def weight_sil(sm1,sm2,sm3,pre_lable1, pre_lable2, pre_lable3):
    sil1 = sil_score(sm1, pre_lable1)
    sil2 = sil_score(sm2, pre_lable2)
    sil3 = sil_score(sm3, pre_lable3)
    W_sum = sil1+sil2+sil3
    w1 = sil1/W_sum
    w2 = sil2/W_sum
    w3 = sil3/W_sum
    return(np.array([w1,w2,w3]))

def weight_PNMI(pre_lable1, pre_lable2, pre_lable3):
    nmi12= normalized_mutual_info_score(pre_lable1, pre_lable2)
    nmi13= normalized_mutual_info_score(pre_lable1, pre_lable3)
    nmi23= normalized_mutual_info_score(pre_lable2, pre_lable3)     
    pnmi1 = 1/2 * (nmi12+nmi13)
    pnmi2 = 1/2 * (nmi12+nmi23)
    pnmi3 = 1/2 * (nmi13+nmi23)
    W_sum = pnmi1+pnmi2+pnmi3
    w1 = pnmi1/W_sum
    w2 = pnmi2/W_sum
    w3 = pnmi3/W_sum
    return(np.array([w1,w2,w3]))

def find_kcluster(k_min,k_max,sm1,sm2,sm3):
    sil_list = []
    for k in range(k_min,k_max+1):
        out1 = sc_pre(sm1,k)
        out2 = sc_pre(sm2,k)
        out3 = sc_pre(sm3,k)
        sil1 = sil_score(sm1,out1)
        sil2 = sil_score(sm2,out2)
        sil3 = sil_score(sm3,out3)
        sil_sum = sil1+sil2+sil3
        sil_list.append(np.array[k,sil_sum])    
    return(sil_list)


def merge_clus(pre_lable1, pre_lable2, pre_lable3,weight_sil,weight_PNMI):
    df_lb = np.vstack((pv_lb,cos_lb))
    df_lb = np.vstack((df_lb,dua_lb))
    df_lb = df_lb.T
    df_lb = np.vstack((df_lb,weight_sil))
    df_lb = np.vstack((df_lb,weight_PNMI))
    df_lb = pd.DataFrame(df_lb)
    return(df_lb)
    