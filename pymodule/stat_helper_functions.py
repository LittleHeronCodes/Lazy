 
from sklearn.metrics import roc_curve, auc

def getRocAUC(test, pred):
	fpr, tpr, _ = roc_curve(test, pvDT.pv_t)
	roc_auc = auc(fpr, tpr)
	return roc_auc

##USAGE

# library(reticulate)
# use_python("/usr/bin/python3.6", required = T)
# # py_config()

# source_python("roc_AUC_function.py")



