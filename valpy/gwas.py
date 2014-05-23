"""
Performs functions necessary for GWAS analysis 
"""


from performetrics import *


def gwas(betaColumn, betaTrueFalse):
	return rmse(betaColumn, betaTrueFalse), mae(betaColumn, betaTrueFalse)