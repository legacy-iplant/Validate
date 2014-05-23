"""
Performs functions necessary for GWAS analysis 
"""


from performetrics import *


def gwas(betaColumn, betaTrueFalse, snpTrueFalse, scoreColumn, threshold):
	return rmse(betaColumn, betaTrueFalse), mae(betaColumn, betaTrueFalse), auc(snpTrueFalse, scoreColumn), tp(snpTrueFalse, threshold, scoreColumn), fp(snpTrueFalse, threshold, scoreColumn), tn(snpTrueFalse, threshold, scoreColumn), fn(snpTrueFalse, threshold, scoreColumn),tpr(snpTrueFalse, threshold, scoreColumn), fpr(snpTrueFalse, threshold, scoreColumn), error(snpTrueFalse, threshold, scoreColumn), sens(snpTrueFalse, threshold, scoreColumn), spec(snpTrueFalse, threshold, scoreColumn), precision(snpTrueFalse, threshold, scoreColumn), youden(snpTrueFalse, threshold, scoreColumn)