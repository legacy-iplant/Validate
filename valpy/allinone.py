""" 
Functions for indentifying and using the command-line to execute Validate for Python
"""


"""Dependencies"""
import getopt, sys
import numpy as np
import pandas as pd
from scipy import stats
import random
import os, data, csv

"""Functions to be used later in the software"""
# Prints introduction graphics for every time the software is run
def initializeGraphics():
	print "###################################################################"
	print "###                                                            ####"
	print "###      Validate for Python!                                  ####"
	print "###      By Dustin A. Landers                                  ####"
	print "###      Contact: (770) 289-8830 -- dustin.landers@gmail.com   ####"
	print "###                                                            ####"
	print "###################################################################"


# Prints all possible command-line arguments to the screen; also ends the execution of the software
def usage():
	print "\n\n\n"
	print "Command-line usage help menu.\n"
	print "--verbose or -v for verbose mode"
	print "--analysis or -a to specify either 'GWAS' or 'prediction' (if blank, Validate assumes GWAS)"
	print "--Folder or -F to input folder of box results (required)"
	print "--Class or -C to specify the known-truth file for used simulation (required)"
	print "--Snp or -S to specify a string for the name of the SNP column in results file (required)"
	print "--Score or -P to specify a string for the name of the scoring column in results file (e.g., p-value; required)"
	print "--beta or -b to specify a string for the name of the estimated SNP effect column in results file"
	print "--severity or -y to specify a severity ratio to use in calculating the H-measure (recommended 1 or pi1/pi0)"
	print "--filename or -f to specify the desired filename for the Validate output file"
	print "--threshold ir -t to specify a desired threshold for classification performetrics where necessary"
	print "--seper or -s to specify either whitespace or comma"
	print "--kttype or -k to specify the type of known-truth file for --class (either OTE or FGS)"
	print "--kttypeseper or -r to specify delimination in known-truth file"
	print "--help or -h to see help menu\n\n"


# Checks for arguments at beginning of the execution of the main function
def checkArgs():
	try:
		opts, args = getopt.getopt(sys.argv[1:], shortopts="vha:F:C:S:P:b:y:f:t:s:k:r", longopts=["verbose", "help", 
			"analysis=", "Folder=", "Class=", "Snp=", "Score=", "beta=", "filename=", "threshold=", "seper=", "kttype=",
			"kttypeseper=", "severity="])

	except getopt.GetoptError as err:
		print(err)
		usage()
		sys.exit()

	# Specifiying initial values of needed variables; unneeded specification when desiring defaults
	verbose = False
	analysis = "GWAS"
	filename = "Results"
	threshold = 0.05
	seper = "whitespace"
	kttype = "OTE"
	kttypeseper = "whitespace"

	# Looping through command-line arguments to replace and/or create initialized values
	for o in opts:
		if o[0] in ("--help", "-h"):
			usage()
			sys.exit()
	for o in opts:
		if o[0] in ("--verbose", "-v"):
			verbose = True
			print ("Verbose mode\n")
	for o in opts:
		if o[0] in ("--Folder", "-F"):
			folder = str(o[1])
			if verbose:
				print "Folder of results files for validation is located in", folder
		if o[0] in ("--analysis", "-a"):
			analysis = str(o[1])
			if verbose:
				print "Analysis method being validated is specified as", analysis
		if o[0] in ("--Class", "-C"):
			truth = str(o[1])
			if verbose:
				print "Truth file is", truth
		if o[0] in ("--Snp", "-S"):
			snp = str(o[1])
			if verbose:
				print "SNP column name in results files is specified as", snp
		if o[0] in ("--Score", "-P"):
			score = str(o[1])
			if verbose:
				print "Scoring column name (e.g., p-value column) in results files is specified as", score
		if o[0] in ("--beta", "-b"):
			beta = str(o[1])
			if verbose:
				print "Estimated SNP Weight column name (e.g., regression betas) in results files is specified as", beta
		if o[0] in ("--filename", "-f"):
			filename = str(o[1])
			if verbose:
				print "Filename specified as", filename
		if o[0] in ("--threshold", "-t"):
			threshold = float(o[1])
			if verbose:
				print "Theshold is set at", threshold
		if o[0] in ("--seper", "-s"):
			seper = str(o[1])
			if verbose:
				print "Delimination of results files is set as", seper
		if o[0] in ("--kttype", "-k"):
			kttype = str(o[1])
			if verbose:
				print "Known-truth data format is set as", kttype
		if o[1] in ("--kttypeseper", "-r"):
			kttypeseper = str(o[1])
			if verbose:
				print "Known-truth data format delimination is set as", kttypeseper
		if o[1] in ("--severity", "-y"):
			severity = float(o[1])
			if verbose:
				print "Severity ratio is specified at", severity

	# Check to see if needed variables are defined
	try:
		folder
	except NameError:
		print "ERROR: Folder of results files to be validated must be specificed."
		usage()
		sys.exit()
	try:
		truth
	except NameError:
		print "ERROR: Known-truth data file must be supplied in order for results to be validated."
		usage()
		sys.exit()
	try:
		snp
	except NameError:
		print "ERROR: Name of SNP column in results files must be specified."
		usage()
		sys.exit()
	try:
		score
	except NameError:
		print "ERROR: Name of scoring column must be specified in order to validate SNP classifications."
		usage()
		sys.exit()

	# Setting beta equal to null if not used; this will placehold the need to not run beta analyses
	try:
		beta
	except NameError:
		beta = None

	# Setting severity equal to null if not used
	try:
		severity
	except NameError:
		severity = None

	return folder, analysis, truth, snp, score, beta, filename, threshold, seper, kttype, kttypeseper, severity

"""
Performs functions necessary for GWAS analysis 
"""


def gwasWithBeta(betaColumn, betaTrueFalse, snpTrueFalse, scoreColumn, threshold):
	return ["h", "rmse", "mae", "r", "r2", "auc", "tp", "fp", "tn", "fn", "tpr", "fpr", "error", "sens", "spec", "precision", "youden"], [h(snpTrueFalse, scoreColumn), rmse(betaColumn, betaTrueFalse), mae(betaColumn, betaTrueFalse), r(betaColumn, betaTrueFalse), r2(betaColumn, betaTrueFalse),auc(snpTrueFalse, scoreColumn), tp(snpTrueFalse, threshold, scoreColumn), fp(snpTrueFalse, threshold, scoreColumn), tn(snpTrueFalse, threshold, scoreColumn), fn(snpTrueFalse, threshold, scoreColumn), tpr(snpTrueFalse, threshold, scoreColumn), fpr(snpTrueFalse, threshold, scoreColumn), error(snpTrueFalse, threshold, scoreColumn), sens(snpTrueFalse, threshold, scoreColumn), spec(snpTrueFalse, threshold, scoreColumn), precision(snpTrueFalse, threshold, scoreColumn), youden(snpTrueFalse, threshold, scoreColumn)]

def gwasWithoutBeta(snpTrueFalse, scoreColumn, threshold):
	return ["h", "auc", "tp", "fp", "tn", "fn", "tpr", "fpr", "error", "sens", "spec", "precision", "youden"], [h(snpTrueFalse, scoreColumn), auc(snpTrueFalse, scoreColumn), tp(snpTrueFalse, threshold, scoreColumn), fp(snpTrueFalse, threshold, scoreColumn), tn(snpTrueFalse, threshold, scoreColumn), fn(snpTrueFalse, threshold, scoreColumn),tpr(snpTrueFalse, threshold, scoreColumn), fpr(snpTrueFalse, threshold, scoreColumn), error(snpTrueFalse, threshold, scoreColumn), sens(snpTrueFalse, threshold, scoreColumn), spec(snpTrueFalse, threshold, scoreColumn), precision(snpTrueFalse, threshold, scoreColumn), youden(snpTrueFalse, threshold, scoreColumn)]

"""
Functions to import both class and results folder files
"""

def getList(folder):
	return os.listdir(folder)


def loadFile(folder, thisFile, seper):
	return data.Data(folder + "/" + thisFile, seper, skiprow=False)


def loadKT(thisFile, seper):
	return data.Data(thisFile, seper, skiprow=True)


def trueFalse(currentSnp, ktSnps):
	if currentSnp in ktSnps:
		return True
	else:
		return False


def writeCSV(filename, keepToWrite, method="wb", exportDelimiter=","):
	with open(filename + ".txt", method) as openFile:
		openFileWriter = csv.writer(openFile, delimiter=exportDelimiter)
		if method == "wb":
			openFileWriter.writerow(keepToWrite[0])
		currentRow = list()
		for item in keepToWrite[1]:
			currentRow.append(item)
		openFileWriter.writerow(currentRow)

"""
Performance measures for testing applications in Validate
"""


def h(snpTrueFalse, scoreColumn):
	n = float(len(scoreColumn))
	n1 = float(sum(snpTrueFalse))
	n0 = n - n1
	pi0 = n0/n
	pi1 = n1/n
	severityRatio = pi1/pi0
	zord = pd.Series(np.ravel(scoreColumn).argsort())
	sc = list()
	count = 0
	for each in scoreColumn:
		sc.append(scoreColumn[zord[count]])
		count += 1
	sc = pd.Series(sc)


def rmse(betaColumn, betaTrueFalse):
	betaColumn = np.array(betaColumn)
	betaTrueFalse = np.array(betaTrueFalse)
	return np.mean(np.square(np.subtract(betaColumn, betaTrueFalse)))


def mae(betaColumn, betaTrueFalse):
	betaColumn = np.array(betaColumn)
	betaTrueFalse = np.array(betaTrueFalse)
	return np.mean(np.absolute(np.subtract(betaColumn, betaTrueFalse)))


def r(betaColumn, betaTrueFalse):
	betaColumn = np.array(betaColumn)
	betaTrueFalse = np.array(betaTrueFalse)
	return stats.stats.pearsonr(betaColumn, betaTrueFalse)[0]


def r2(betaColumn, betaTrueFalse):
	betaColumn = np.array(betaColumn)
	betaTrueFalse = np.array(betaTrueFalse)
	return np.square(stats.stats.pearsonr(betaColumn, betaTrueFalse)[0])


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


def tp(snpTrueFalse, threshold, scoreColumn):
	testColumn = list()
	for each in scoreColumn:
		if float(each) < threshold:
			testColumn.append(True)
		else:
			testColumn.append(False)
	count = 0
	truePositives = 0
	for each in testColumn:
		if each is True and snpTrueFalse[count] is True:
			truePositives += 1
		count += 1
	return truePositives


def fp(snpTrueFalse, threshold, scoreColumn):
	testColumn = list()
	for each in scoreColumn:
		if float(each) < threshold:
			testColumn.append(True)
		else:
			testColumn.append(False)
	count = 0
	falsePositives = 0
	for each in testColumn:
		if each is True and snpTrueFalse[count] is False:
			falsePositives += 1
		count += 1
	return falsePositives


def tn(snpTrueFalse, threshold, scoreColumn):
	testColumn = list()
	for each in scoreColumn:
		if float(each) < threshold:
			testColumn.append(True)
		else:
			testColumn.append(False)
	count = 0
	trueNegatives = 0
	for each in testColumn:
		if each is False and snpTrueFalse[count] is False:
			trueNegatives += 1
		count += 1
	return trueNegatives


def fn(snpTrueFalse, threshold, scoreColumn):
	testColumn = list()
	for each in scoreColumn:
		if float(each) < threshold:
			testColumn.append(True)
		else:
			testColumn.append(False)
	count = 0
	falseNegatives = 0
	for each in testColumn:
		if each is False and snpTrueFalse[count] is True:
			falseNegatives += 1
		count += 1
	return falseNegatives


def tpr(snpTrueFalse, threshold, scoreColumn):
	truePositives = tp(snpTrueFalse, threshold, scoreColumn)
	count = 0.0
	for each in snpTrueFalse:
		if each is True:
			count += 1.0
	return float(truePositives/count)


def fpr(snpTrueFalse, threshold, scoreColumn):
	falsePositives = fp(snpTrueFalse, threshold, scoreColumn)
	count = 0.0
	for each in snpTrueFalse:
		if each is False:
			count += 1.0
	return float(falsePositives/count)


def error(snpTrueFalse, threshold, scoreColumn):
	truePositives = float(tp(snpTrueFalse, threshold, scoreColumn))
	falsePositives = float(fp(snpTrueFalse, threshold, scoreColumn))
	trueNegatives = float(tn(snpTrueFalse, threshold, scoreColumn))
	falseNegatives = float(fn(snpTrueFalse, threshold, scoreColumn))
	return (falseNegatives + falsePositives) / (truePositives + trueNegatives + falsePositives + falseNegatives)


def sens(snpTrueFalse, threshold, scoreColumn):
	truePositives = float(tp(snpTrueFalse, threshold, scoreColumn))
	falseNegatives = float(fn(snpTrueFalse, threshold, scoreColumn))
	return truePositives / (truePositives + falseNegatives)


def spec(snpTrueFalse, threshold, scoreColumn):
	trueNegatives = float(tn(snpTrueFalse, threshold, scoreColumn))
	falsePositives = float(fp(snpTrueFalse, threshold, scoreColumn))
	return trueNegatives / (trueNegatives + falsePositives)


def precision(snpTrueFalse, threshold, scoreColumn):
	truePositives = float(tp(snpTrueFalse, threshold, scoreColumn))
	falsePositives = float(fp(snpTrueFalse, threshold, scoreColumn))
	return truePositives / (truePositives + falsePositives)


def youden(snpTrueFalse, threshold, scoreColumn):
	sensitivity = float(sens(snpTrueFalse, threshold, scoreColumn))
	specificity = float(spec(snpTrueFalse, threshold, scoreColumn))
	return sensitivity + sensitivity - 1.0


"""Main function and execution"""
def main():
	initializeGraphics()
	folder, analysis, truth, snp, score, beta, filename, threshold, seper, kttype, kttypeseper, severity = checkArgs()
	appOutputList = checkList(getList(folder))
	ktFile = loadKT(truth, kttypeseper)


	if kttype == "OTE":
		acquiredData = loadFile(folder, appOutputList[0], seper)
		snpColumnNo = acquiredData.header.index(snp)
		snpColumn = list()
		for each in acquiredData.data.iteritems():
			snpColumn.append(each[1][snpColumnNo])
		
		ktSnps = list()
		for each in ktFile.data.iteritems():
			ktSnps.append(each[1][0])
		ktBetas = list()
		for each in ktFile.data.iteritems():
			ktBetas.append(each[1][1])

		snpTrueFalse = list()
		for each in snpColumn:
			snpTrueFalse.append(trueFalse(each, ktSnps))
		
		if beta is not None:
			betaTrueFalse = list()
			count = 0
			for each in snpTrueFalse:
				if each is True:
					current = snpColumn[count]
					match = ktSnps.index(current)
					thisBeta = ktBetas[match]
					betaTrueFalse.append(float(thisBeta))
				else:
					betaTrueFalse.append(float(0))
				count += 1

		if severity is None:
			severity = float(len(ktSnps))/float(len(snpTrueFalse) - len(ktSnps))

	firstForHeader = True
	for each in appOutputList:
		acquiredData = loadFile(folder, each, seper)
		snpColumnNo = acquiredData.header.index(snp)
		snpColumn = list()
		for each in acquiredData.data.iteritems():
			snpColumn.append(each[1][snpColumnNo])

		scoreColumnNo = acquiredData.header.index(score)
		scoreColumn = list()
		for each in acquiredData.data.iteritems():
			scoreColumn.append(float(each[1][scoreColumnNo]))

		if beta is not None:
			betaColumnNo = acquiredData.header.index(beta)
			betaColumn = list()
			for each in acquiredData.data.iteritems():
				betaColumn.append(float(each[1][betaColumnNo]))

		if analysis == "GWAS" and firstForHeader:
			if beta is not None:
				keepToWrite = gwasWithBeta(betaColumn, betaTrueFalse, snpTrueFalse, scoreColumn, threshold)
				writeCSV(filename, keepToWrite, "wb", "\t")
			if beta is None:
				keepToWrite = gwasWithoutBeta(snpTrueFalse, scoreColumn, threshold)
				writeCSV(filename, keepToWrite, "wb", "\t")
		else:
			if beta is not None:
				keepToWrite = gwasWithBeta(betaColumn, betaTrueFalse, snpTrueFalse, scoreColumn, threshold)
				writeCSV(filename, keepToWrite, "a", "\t")
			if beta is None:
				keepToWrite = gwasWithoutBeta(snpTrueFalse, scoreColumn, threshold)
				writeCSV(filename, keepToWrite, "a", "\t")
		firstForHeader = False


if __name__ == "__main__":
	main()