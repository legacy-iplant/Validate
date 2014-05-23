"""
Performs functions necessary for GWAS analysis 
"""


from performetrics import *


def gwas(betaColumn, betaTrueFalse, snpTrueFalse, scoreColumn):
	return rmse(betaColumn, betaTrueFalse), mae(betaColumn, betaTrueFalse), auc(snpTrueFalse, scoreColumn)