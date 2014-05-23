"""
Performs functions necessary for GWAS analysis 
"""


from performetrics import *


def gwas(betaColumn, betaTrueFalse, snpTrueFalse, scoreColumn, threshold):
	return rmse(betaColumn, betaTrueFalse), mae(betaColumn, betaTrueFalse), auc(snpTrueFalse, scoreColumn), tp(snpTrueFalse, threshold, scoreColumn), fp(snpTrueFalse, threshold, scoreColumn)