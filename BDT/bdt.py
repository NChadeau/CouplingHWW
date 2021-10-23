import Channels
import pandas as pd
import numpy as np
import math
import ROOT
import sys
import copy
import argparse
import psutil
import matplotlib.pyplot as plt
import matplotlib as mpl
from joblib import dump, load

from sklearn.tree import export_graphviz
from sklearn.metrics import classification_report, roc_auc_score
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

inputNames = ["reco_E4jets","reco_Pt4jets","reco_M4jets","reco_WBigMass_Energy","reco_WBigMass_Pt","reco_WBigMass_Mass","reco_WBigMass_CosJets","reco_WSmallMass_Energy","reco_WSmallMass_Pt","reco_WSmallMass_Mass","reco_WSmallMass_CosJets","reco_Cos","reco_njets","reco_Y12","reco_Y13","reco_Y14","reco_Y23","reco_Y24","reco_Y34"]

#inputNames = ['reco_E4jets','reco_Pt4jets', 'reco_M4jets', 'reco_njets', 'reco_Y12']

signalBranches = ["mc_HiggsDecay", "mc_DecayProcess"]
cuts = ["reco_BoolDumpEvent"]

def readFile(runNumber, isSignal, lumi, ePol, pPol) :
	print(f'Read file {runNumber}...')
	root_file = ROOT.TFile(f'/home/nicolas/Documents/ep->nnH_(H->WW->qqqq)/MyAnalysis/script/{runNumber}.root')
	tree = root_file.Get("tree")
	
	df = ROOT.RDataFrame(tree)
	X = df.AsNumpy(columns=inputNames + cuts)
	
	data = pd.DataFrame()
	for name in inputNames + cuts:
		data[name] = X[name]
	
	nEventsInFile = data.shape[0]
	
	is_NaN = data.isnull()
	row_has_NaN = is_NaN.any(axis=1)

	eventsToRemove = row_has_NaN.values

	targetArray = np.full(nEventsInFile, 0)

	if isSignal : 
		temp = df.AsNumpy(columns=signalBranches)
		arrayOfSignal = ((temp["mc_HiggsDecay"] == 2424) & (temp["mc_DecayProcess"] == 808))
		
		targetArray[arrayOfSignal] = 1

	data['target'] = targetArray
	data['target'] = data['target'].astype('category')
	
	#processIdArray = np.full(nEventsInFile, str(runNumber))
	data['processID'] = str(runNumber)
	data['processID'] = data['processID'].astype('category')
	
	weight = Channels.getWeight(runNumber, nEventsInFile, lumi, ePol, pPol)
	#data['weight'] = np.full(nEventsInFile, weight)
	data['weight'] = weight
	data['weight'] = data['weight'].astype('float32')
	
	eventsToRemove = eventsToRemove | (data['reco_BoolDumpEvent'] == 1)
	data = data[~eventsToRemove]
	
	root_file.Close()
	
	data = data.drop(columns=cuts)
	
	print(f'{nEventsInFile} events in file, {data.shape[0]} events kept, {data.shape[0]*weight:.2f} in reality')
	if isSignal :
		print(f'{np.sum(targetArray)} signal events, {np.sum(targetArray)*weight:.2f} in reality')
	
	return data


def compare_train_test(scores_train, scores_test, y_train, y_test, bins=100):

	decisions = []
	for score, y in ((scores_train, y_train), (scores_test, y_test)):
		d1 = score[y > 0.5]
		d2 = score[y < 0.5]
		decisions += [d1, d2]

	low = min(np.quantile(d, 0.05) for d in decisions)
	high = max(np.quantile(d, 0.999) for d in decisions)

	low_high = (low, high)

	plt.hist(decisions[0], color='r', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True, label='S (train)')
	plt.hist(decisions[1], color='b', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True, label='B (train)')

	hist, bins = np.histogram(decisions[2], bins=bins, range=low_high, density=True)
	scale = len(decisions[2]) / sum(hist)
	err = np.sqrt(hist * scale) / scale

	width = (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	plt.errorbar(center, hist, yerr=err, fmt='o', c='r', label='S (test)')

	hist, bins = np.histogram(decisions[3], bins=bins, range=low_high, density=True)
	scale = len(decisions[2]) / sum(hist)
	err = np.sqrt(hist * scale) / scale

	plt.errorbar(center, hist, yerr=err, fmt='o', c='b', label='B (test)')

	plt.xlabel("BDT output")
	plt.ylabel("Arbitrary units")
	plt.legend(loc='best')

	plt.show()


if __name__ == '__main__' :

	parser = argparse.ArgumentParser()
	parser.add_argument('-d','--data', help='Data file', required=False)
	parser.add_argument('-e','--ePol', help='E polarisation', required=False)
	parser.add_argument('-p','--pPol', help='P polarisation', required=False)
	parser.add_argument('-l','--lumi', help='Integrated Luminosity', required=False)
	parser.add_argument('-m','--model', help='Load existing estimator', required=False)
	args = vars(parser.parse_args())

	loadFile = False
	if args['data'] :
		loadFile = args['data']

	ePol = -0.8
	pPol = 0.3
	lumi = 900
	
	if args['ePol'] :
		ePol = float(args['ePol'])
	if args['pPol'] :
		pPol = float(args['pPol'])
	if args['lumi'] :
		lumi = float(args['lumi'])
	
	print(f'ePol : {ePol}, pPol : {pPol}, lumi : {lumi}')
	
	signal = ['402007','402008']
	higgs = ['402001','402002','402003','402004','402005','402006','402009','402010','402011','402012']
	fermions2 = ['500006','500008','500010','500012']
	fermions4 = ['500070','500072','500074','500076','500082','500084','500098','500100','500101','500102','500103','500104','500105','500106','500107','500108','500110','500112','500113','500114','500115','500116','500117','500118','500119','500120','500122','500124']
	
	dataFrame = pd.DataFrame()
	# read files
	if not loadFile:
		print('Data file not provided : loading from ROOT files...')
		dataFrameList = []
		for processID, osef in Channels.channels.items() :
			isSignal = False
			if processID in signal :
				isSignal = True
				
			data = readFile(processID, isSignal, lumi, ePol, pPol)
			dataFrameList.append(data)
		
		dataFrame = pd.concat(dataFrameList)
		dataFrameList.clear()
		dataFrame.to_pickle(f'e{ePol}_p{pPol}_{lumi}.pkl')
	else:
		print(f'Loading events from {loadFile}...')
		dataFrame = pd.read_pickle(loadFile)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	cax = ax.matshow(dataFrame.loc[:,inputNames].corr(), interpolation='nearest')
	fig.colorbar(cax)


	xaxis = np.arange(len(inputNames))
	ax.set_yticks(xaxis)
	ax.set_yticklabels(inputNames)
	
	plt.xticks([])
	plt.show()

	signalTestProp = 0.5
	bkgTestNumber = 3 #times number of signal events
	
	#y = dataFrame['target']
	#data_train, data_test = train_test_split(dataFrame, train_size=signalTestProp, random_state=42, stratify=y)
	data_train, data_test = train_test_split(dataFrame, train_size=signalTestProp, random_state=42)

	print('compute weights for train...')
	signalTotalWeight = data_train.loc[data_train['target'] == 1, 'weight'].sum()
	bkgTotalWeight = data_train.loc[data_train['target'] == 0, 'weight'].sum()
	cheatWeight = bkgTotalWeight/signalTotalWeight
	data_train.loc[data_train['target'] == 1 , 'weight'] *= cheatWeight

	osef = [dataFrame]
	del dataFrame
	del osef
	print('original dataFrame released...')
	
	bdt = None

	if args['model'] :
		bdt = load(args['model'])
	else :
		#dt = DecisionTreeClassifier(max_depth=2)
		#bdt = AdaBoostClassifier(algorithm='SAMME.R', n_estimators=1000, learning_rate=0.5)
		bdt = GradientBoostingClassifier(n_estimators=500, learning_rate=1.0, max_depth=2, random_state=0, verbose=1)
		#bdt = RandomForestClassifier(max_depth=1, random_state=0, n_estimators=208, n_jobs=-1, verbose=2)
	
		bdt.fit(data_train[inputNames], data_train['target'], sample_weight=data_train['weight'])
		dump(bdt, 'filename.joblib')

	importances = bdt.feature_importances_
	for i,v in enumerate(importances):
		print(f'{inputNames[i]} : {v*100:.6f} %')
		
	countPerProcessID = data_test.groupby(['processID']).size()
	for processID, nEvents in countPerProcessID.items() :
		weight = Channels.getWeight(processID, nEvents, lumi, ePol, pPol)
		data_test.loc[data_test['processID'] == processID, 'weight'] = weight
		print(f'{processID} : {data_test.loc[data_test["processID"] == processID, "weight"].mean():.2f}')
		print(f'{processID} : {nEvents} events in test, {nEvents*weight:.2f} in reality')
		
	#print(f'Total events : {data_test.groupby(["target"])["weight"].sum()}')
	#print(f'Total events : {data_test.groupby(["target","processID"]).size()}')
	
	print('compute scores...')
	y_pred = bdt.predict(data_test[inputNames])
	
	data_test.to_pickle('/tmp/data_test')
	osef = [data_test]
	del data_test
	del osef
	
	scores_train = bdt.decision_function(data_train[inputNames])
	
	y_train = data_train['target'].values
	data_train.to_pickle('/tmp/data_train')
	osef = [data_train]
	del data_train
	del osef
	
	data_test = pd.read_pickle('/tmp/data_test')
	scores_test = bdt.decision_function(data_test[inputNames])
	
	y_test = data_test['target'].values
	#scores_train = bdt.predict_proba(X_train)[:,1]
	#scores_test = bdt.predict_proba(X_test)[:,1]

	print('scores computed...')
	
	def computeSignificance(scoreCut, scores, realClasses, weights) :
		toSelect = scores > scoreCut
		selected = realClasses[toSelect]
		weightsSelected = weights[toSelect]
		
		totalSelected = np.sum(weightsSelected)
		signalSelected = np.sum(selected * weightsSelected)

		significance = signalSelected/math.sqrt(totalSelected)
		return significance
		
	minScore = np.quantile(scores_test, 0.05)
	maxScore = np.quantile(scores_test, 0.999)
	
	scoresToTest = np.linspace(minScore, maxScore, 100, endpoint=False)
	significances = []
	for scoreCut in scoresToTest :
		significance = computeSignificance(scoreCut, scores_test, data_test['target'], data_test['weight'])
		significances.append(significance)
	
	positionOfMaxSignificance = np.argmax(significances)
	bestCut = scoresToTest[positionOfMaxSignificance]
	y_bestCut = 1*(scores_test > bestCut)

	plt.plot(scoresToTest, significances, 'bo')
	plt.xlabel("BDT output")
	plt.ylabel("Significance")
	plt.show()

	print(classification_report(data_test['target'], y_bestCut, target_names=['bkg','signal'], sample_weight = data_test['weight']))
	#export_graphviz(dt, out_file="test.dot", feature_names=inputNames, class_names=['bkg','signal'], rounded=True, filled=True)
	
	def getChannelType(processID, target):
		
		channelType = np.full(processID.shape[0], 'unknown')
		channelType[target==1] = 'signal'
		
		pInHiggs1 = np.array([(lambda p: p in signal)(p) for p in processID])
		pInHiggs1 = (pInHiggs1) & (target==0)
		pInHiggs2 = [(lambda p: p in higgs)(p) for p in processID]
		pInFermions2 = [(lambda p: p in fermions2)(p) for p in processID]
		pInFermions4 = [(lambda p: p in fermions4)(p) for p in processID]
		
		channelType[pInHiggs1] = 'Higgs 1'
		channelType[pInHiggs2] = 'Higgs 2'
		channelType[pInFermions2] = '2 fermi'
		channelType[pInFermions4] = '4 fermi'

		return channelType

	data_test['pred'] = y_bestCut
	data_test['channelType'] = getChannelType(data_test['processID'], data_test['target'])
	
	print(data_test['channelType'])
	
	finalSelection = data_test.groupby(['channelType','pred'])['weight']
	#print(finalSelection.groups.keys())
	
	selectedSignal = 0
	selectedAll = 0
	for channelType in ['signal', 'Higgs 1', 'Higgs 2', '2 fermi', '4 fermi']:
		selected = finalSelection.get_group((channelType, 1)).sum()
		if channelType == 'signal' :
			selectedSignal = selected
		selectedAll += selected
		rejected = finalSelection.get_group((channelType, 0)).sum()
		finalEfficiency = selected/(selected+rejected)
		print(f'{channelType} cross section : {(selected+rejected)/lumi:.6f}, {channelType} expected : {selected+rejected:.6f} , {channelType} selected : {selected:.6f}')
		print(f'{channelType} : {100*finalEfficiency:.6f} %')
	maxSignificance = np.max(significances)
	number = np.where(significances == maxSignificance)
	BDToutput = scoresToTest[number[0][0]]
	purity = selectedSignal/selectedAll
	print(f'BDT output : {BDToutput:.6f}')
	print(f'Maximum significance : {maxSignificance:.6f}')
	print(f'Purity : {100*purity:.6f} %')
	print(f'Relative uncertainty on the cross section: {100.0/maxSignificance:.6f} %')

	compare_train_test(scores_train, scores_test, y_train, y_test)
