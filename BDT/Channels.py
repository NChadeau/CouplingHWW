channels = {
    '402007': {'xSect': 60.3509, 'pol': 'LR'},
    '402008': {'xSect': 21.4642, 'pol': 'RL'},
    '402001': {'xSect': 17.6715, 'pol': 'LR'},
    '402002': {'xSect': 11.1389, 'pol': 'RL'},
    '402003': {'xSect': 16.9707, 'pol': 'LR'},
    '402004': {'xSect': 10.8691, 'pol': 'RL'},
    '402005': {'xSect': 16.9407, 'pol': 'LR'},
    '402006': {'xSect': 10.8434, 'pol': 'RL'},
    '402009': {'xSect': 67.1107, 'pol': 'LR'},
    '402010': {'xSect': 42.9279, 'pol': 'RL'},
    '402011': {'xSect': 343.030, 'pol': 'LR'},
    '402012': {'xSect': 219.486, 'pol': 'RL'},
    '500006': {'xSect': 21214, 'pol': 'LR'},
    '500008': {'xSect': 16363, 'pol': 'RL'},
    '500010': {'xSect': 127966, 'pol': 'LR'},
    '500012': {'xSect': 70416.7, 'pol': 'RL'},
    '500070': {'xSect': 12389.3, 'pol': 'LR'},
    '500072': {'xSect': 225.569, 'pol': 'RL'},
    '500074': {'xSect': 838.079, 'pol': 'LR'},
    '500076': {'xSect': 466.816, 'pol': 'RL'},
    '500082': {'xSect': 18779.1, 'pol': 'LR'},
    '500084': {'xSect': 173.468, 'pol': 'RL'},
    '500098': {'xSect': 1637.06, 'pol': 'LR'},
    '500100': {'xSect': 55.4225, 'pol': 'RL'},
    '500101': {'xSect': 1155.83, 'pol': 'LL'},
    '500102': {'xSect': 1423.31, 'pol': 'LR'},
    '500103': {'xSect': 1157.2, 'pol': 'RR'},
    '500104': {'xSect': 1219.4, 'pol': 'RL'},
    '500105': {'xSect': 190.531, 'pol': 'LL'},
    '500106': {'xSect': 10264, 'pol': 'LR'},
    '500107': {'xSect': 190.637, 'pol': 'RR'},
    '500108': {'xSect': 86.6962, 'pol': 'RL'},
    '500110': {'xSect': 453.87, 'pol': 'LR'},
    '500112': {'xSect': 131.22, 'pol': 'RL'},
    '500113': {'xSect': 5619, 'pol': 'LL'},
    '500114': {'xSect': 5774.74, 'pol': 'LR'},
    '500115': {'xSect': 5629.24, 'pol': 'RR'},
    '500116': {'xSect': 5675.39, 'pol': 'RL'},
    '500117': {'xSect': 63.521, 'pol': 'LL'},
    '500118': {'xSect': 3421.97, 'pol': 'LR'},
    '500119': {'xSect': 63.4832, 'pol': 'RR'},
    '500120': {'xSect': 29.4883, 'pol': 'RL'},
    '500122': {'xSect': 195.126, 'pol': 'LR'},
    '500124': {'xSect': 40.683, 'pol': 'RL'}
}

def getPolarisationWeights(ePol, pPol):

	eR = 0.5*(ePol+1)
	eL = 0.5*(1-ePol)
	pR = 0.5*(pPol+1)
	pL = 0.5*(1-pPol)
	
	weights = {'LR' : eL*pR, 'RL' : eR*pL, 'LL' : eL*pL, 'RR' : eR*pR}
	return weights


def getWeight(processID, nEvents, luminosity, ePol, pPol) :
	nEventsExpected = channels[f'{processID}']['xSect'] * luminosity
	pol = channels[f'{processID}']['pol']
	polWeight = getPolarisationWeights(ePol, pPol)[pol]
	nEventsExpected *= polWeight
	
	return nEventsExpected/nEvents
