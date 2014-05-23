"""
Performance measures for testing Validate
"""


import numpy as np
from scipy import stats


def rmse(betaColumn, betaTrueFalse):
	betaColumn = np.array(betaColumn)
	betaTrueFalse = np.array(betaTrueFalse)
	return np.mean(np.square(np.subtract(betaColumn, betaTrueFalse)))


def mae(betaColumn, betaTrueFalse):
	betaColumn = np.array(betaColumn)
	betaTrueFalse = np.array(betaTrueFalse)
	return np.mean(np.absolute(np.subtract(betaColumn, betaTrueFalse)))


def auc(snpTrueFalse, scoreColumn):
	scoreColumn = np.array(scoreColumn)
	snpTrueFalse = np.array(snpTrueFalse)
	x1 = scoreColumn[snpTrueFalse == True]
	n1 = x1.size
	x2 = scoreColumn[snpTrueFalse == False]
	n2 = x2.size
	r = stats.rankdata(np.hstack((x1,x2)))
	auc = (np.sum(r[0:n1]) - n1 * (n1+1)/2) / (n1 * n2)
	return 1 - auc