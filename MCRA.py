from __future__ import print_function
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
import csv
#import json
from flask import jsonify
import sqlite3
import cPickle

########## Populating arrays from CSV files ########## 

def main():
	smileStr=raw_input()
	pains(smileStr)

def pains(smileStr):
	painsData = []	#tuples of data about everyPAINS alert taken from master CSV

	with open('PAINS_master.csv') as painsAlertCSV:
		reader = csv.DictReader(painsAlertCSV, delimiter=',')
		for row in reader:
			tempDataTuple =(row['SMARTS'], row['\xef\xbb\xbfALERT NAME'], row['NPubChem'], row['NDCM'], row['Luciferase'], row['\xce\xb2-lactamase'], row['Fluorescence'], row['All Assays']) #SMARTS string for each alert is at index 0
			painsData.append(tempDataTuple)
		
	painsCompoundData = [] #array of tuples, everything about compounds, to be used for MACCSkeys analysis
	nonPainsCompoundData = []

	with open('PAINS_compounds_noDuplicates.csv') as painsCompoundsCSV:
		reader = csv.DictReader(painsCompoundsCSV, delimiter=',')
		for row in reader:
			tempDataTuple =((row['\xef\xbb\xbfCID'], row['SMILES'], row['AllAssays_Active'], row['AllAssays_Total'], row['LuciferaseAssays_Active'], row['LuciferaseAssays_Total'], row['BetaLactamaseAssays_Active'], row['BetaLactamaseAssays_Total'], row['FluorescenceAssays_Active'], row['FluorescenceAssays_Total'], row['PAINS_A'], row['PAINS_B'], row['PAINS_C'])) #3 PAINS alerts b/c a compound can have more than 1
			painsCompoundData.append(tempDataTuple)
			
	with open('PAINS_compounds_noDuplicates.csv') as nonPainsCompoundsCSV:
		reader = csv.DictReader(nonPainsCompoundsCSV, delimiter=',')
		for row in reader:
			tempDataTuple =((row['\xef\xbb\xbfCID'], row['SMILES'], row['AllAssays_Active'], row['AllAssays_Total'], row['LuciferaseAssays_Active'], row['LuciferaseAssays_Total'], row['BetaLactamaseAssays_Active'], row['BetaLactamaseAssays_Total'], row['FluorescenceAssays_Active'], row['FluorescenceAssays_Total']))
			nonPainsCompoundData.append(tempDataTuple)

	########## Flagging PAINS Alerts on query compound ##########

	matchIndices = {}	#indices of flagged alert on molecule (for highlighting)
	flaggedAlerts = []	#information about each FLAGGED PAINS alert, taken from painsData[]

	querySmiles = str(smileStr) #SMILE string passed in through POST request.
	queryMol = Chem.MolFromSmiles(querySmiles) #Convert queried SMILE string to a mol.
	if(queryMol==None): return 'Error, poorly formed smile string.'
	for alert in painsData:
		painsAlert = Chem.MolFromSmarts(alert[0]) #converting alert SMARTS string to mol (SMARTS located at index 0 of alert tuple)
		print('HELLOOJIODFJIOJEOIJFIODJOIDFJIO')
		if Chem.AddHs(queryMol).HasSubstructMatch(painsAlert): #search with explicit hydrogens added
			matches = []
			matches.extend(Chem.AddHs(queryMol).GetSubstructMatches(painsAlert)) #Create a dictionary of the name of the alert and the matching substructre indicies.
			matchDic[alert[1]] = matches
			matchIndices.append(matchDic) #adds indices of PAINS alert substructure on query molecule 
			flaggedAlerts.append(alert)
		
		elif queryMol.HasSubstructMatch(painsAlert): #search without explicit hydrogens added
			matches = []
			matches.extend(queryMol.GetSubstructMatches(painsAlert)) #Create a dictionary of the name of the alert and the matching substructre indicies.
			matchDic[alert[1]] = matches
			#matchDic[alert[1]] = queryMol.GetSubstructMatches(painsAlert)
			matchIndices.append(matchDic)
			flaggedAlerts.append(alert)

		
	print('matches for input SMILES')
	for match in flaggedAlerts:
		print(match)
	
	########## Getting Compounds with same alert ##########

	sameAlertCompounds = [] #all compounds with same alert as query molecule

	for painsAlert in flaggedAlerts:
		for row in painsCompoundData: #10,11,12 are where the PAINS alerts are in the tuples of flaggedAlerts
			if row[10] == painsAlert[1]: #painsAlert[1] is the alert name
				sameAlertCompounds.append([row,painsAlert])
			if row[11] == painsAlert[1]: 
				sameAlertCompounds.append([row,painsAlert])
			if row[12] == painsAlert[1]: 
				sameAlertCompounds.append([row,painsAlert])
	
	########## MACCS Fingerprint analysis ########## 

	mostSimilarTc = 0.0		#Tc is a measure of how similar the compound is to the query molecule
	painsMostSimilarDic = None	#compound with highest Tc value
	painsSimilarDics = []	#all compounds with Tc > 0.9, not including the mostSimilarMol

	queryFPS = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(querySmiles)) #FPS=fingerprints which is how MACCSkeys analysis works

	for compound in sameAlertCompounds:
		compoundMol = Chem.MolFromSmiles(compound[0][1])
		painsSubstructure =  Chem.MolFromSmarts(compound[1][0])
		highlights = compoundMol.GetSubstructMatch(painsSubstructure)
		compoundFPS = AllChem.GetMorganFingerprint(compoundMol,2) #1 is index of the SMILES string in the tuple 'compound' in sameAlertCompounds
		compoundTc = DataStructs.FingerprintSimilarity(compoundFPS, queryFPS)
		compoundDic = { 'compound' : compound[0], 'tc' : compoundTc, 'highlights' : highlights }
		if compoundTc > mostSimilarTc : #finding mostSimilarMol
			if(mostSimilarTc > 0.9): painsSimilarMols.append(painsMostSimilarMol) #append older mostSimilar if above 0.9
			mostSimilarTc = compoundTc #replace the old best tc with the new best similarity.
			painsMostSimilarDic = compoundDic
		elif compoundTc > 0.9: #Tc value > 0.9
			#mostSimilarMol = compoundDic
			painsSimilarDics.append(compoundDic)

	mostSimilarTc = 0.0	
	nonPainsMostSimilarMol = None
	nonPainsSimilarMols = []
	

	ls=[]
	ls.append({'most_similar':painsMostSimilarDic})
	ls.append({'other_similar_mols' : painsSimilarDics})
#	ls.append({'non_PAINS_most_similar':nonPainsMostSimilarMol})
#	ls.append({'non_PAINS_other_similar_mols' : nonPainsSimilarMols})
	ls.append({'match_indices' : matchIndices})
	ls.append({'flagged_alerts' : flaggedAlerts})
	
	similar_mols_json = jsonify(ls)
	return similar_mols_json


if __name__ == '__main__':
	main()

	
	
	
	
	
	
	
	
	
	
	
from rdkit.Chem import AllChem	
	
def pains2(smileStr):
	painsData = []
	
	with open('PAINS_master.csv') as painsAlertCSV:
		reader = csv.DictReader(painsAlertCSV, delimiter=',')
		for row in reader:
			tempDataTuple =(row['SMARTS'], row['\xef\xbb\xbfALERT NAME'], row['NPubChem'], row['NDCM'], row['Luciferase'], row['\xce\xb2-lactamase'], row['Fluorescence'], row['All Assays']) #SMARTS string for each alert is at index 0
			painsData.append(tempDataTuple)
	
	########## Flagging PAINS Alerts on query compound ##########

	matchIndices = {}	#indices of flagged alert on molecule (for highlighting)
	flaggedAlerts = []	#information about each FLAGGED PAINS alert, taken from painsData[]

	querySmiles = str(smileStr) #SMILE string passed in through POST request.
	queryMol = Chem.MolFromSmiles(querySmiles) #Convert queried SMILE string to a mol.
	if(queryMol==None): return 'Error, poorly formed smile string.'
	
	for alert in painsData:
		painsAlert = Chem.MolFromSmarts(alert[0]) #converting alert SMARTS string to mol (SMARTS located at index 0 of alert tuple)
	
		if Chem.AddHs(queryMol).HasSubstructMatch(painsAlert): #search with explicit hydrogens added
			all_matches = Chem.AddHs(queryMol).GetSubstructMatches(painsAlert) #adds indices of PAINS alert substructure on query molecule 
			matches = []
			for i in range(len(all_matches)):
				b = True
				for j in range(i):
					diff = 0
					for x, y in zip(all_matches[i], all_matches[j]):
						diff = diff + abs(x - y)
					if(diff > 3): b = False #diff is very large, two highlights are similar
				if(b): matches.append(all_matches[i])
			matchIndices[alert[1]] =  matches
			flaggedAlerts.append(alert)
		elif queryMol.HasSubstructMatch(painsAlert): #search without explicit hydrogens added
			matchIndices[alert[1]] = queryMol.GetSubstructMatches(painsAlert)
			flaggedAlerts.append(alert)

		
	print('PAINS matches for input SMILES')
	wherePains = ''
	orString = ''
	painsSQLQueries = []
	hasPains = True
	for i in range(len(flaggedAlerts)):
		match = flaggedAlerts[i]
		print(match)
		tmp_str = "PAINS_A = '"+ match[1] + "' OR PAINS_B = '"+ str(match[1]) + "' OR PAINS_C = ' " + str(match[1]) + " ' "
		painsSQLQueries.append("SELECT CID, MACCS_fingerprint FROM painsCompounds WHERE " + tmp_str + " ; ")
		wherePains = wherePains + orString + tmp_str
		#builds a global query to be inverted to get all matches without same PAINS as the query.
		orString = ' OR '
	if(len(flaggedAlerts) == 0):
		wherePains = '1<>1'
		hasPains = False
		print("No PAINS")
	########## Getting Compounds with same alert ##########

	sameAlertCompounds = [] #all compounds with same alert as query molecule
	
	#sqliteQuery1 = 
	
	sqliteQuery2 = "SELECT CID, MACCS_fingerprint FROM painsCompounds WHERE NOT (" + wherePains +");"
	
	sqliteQuery3 = "SELECT CID, MACCS_fingerprint FROM nonPainsCompounds;"
	
	conn = sqlite3.connect("pains.db")
	#conn.text_factory = str
	
	cur = conn.cursor()
	cur2 = conn.cursor()
	queryFPS = AllChem.GetMorganFingerprint(queryMol,2)

#	cur.execute(sqliteQuery1)
#	rows = cur.fetchall()	
#	hasPains = True
#	if(len(rows) == 0): #15 from anywhere	
#		print("No PAINS")
#		hasPains = False
#	else: #5 from A, 5 from B, 5 from C, 5 from Other, 5 from Non-pains
#		print("Has PAINS")
	
	#SMILES varchar(255), PAINS_A varchar(31),PAINS_B varchar(31), PAINS_C varchar(31),AllAssays_Active integer, AllAssays_Total integer,LuciferaseAssays_Active integer, 
	#LuciferaseAssays_Total integer, BetaLactamaseAssays_Active integer, BetaLactamaseAssays_Total integer, FluorescenceAssays_Active integer, FluorescenceAssays_Total integer, 
	#PAINS_A_highlights blob, PAINS_B_highlights blob, PAINS_C_highlights 
	
	def rowToDic(cid, tc, row, pains, pains_highlights):
		dic = {'cid': cid, 'tc' : tc, 'smile' : row[0]}
		#'pains_A' : row[1], 'pains_B' : row[2], 'pains_C': row[3], 'pains_A_highlights' : row[12], 'pains_B_highlights' : row[13], 'pains_C_highlights' : row[14]}
		dic['pains'] = pains
		dic['pains_highlights'] = pains_highlights
		dic['AllAssays_Active'] = row[4]
		dic['AllAssays_Total'] = row[5]
		dic['LuciferaseAssays_Active'] = row[6]
		dic['LuciferaseAssays_Total'] = row[7]
		dic['BetaLactamaseAssays_Active'] = row[8]
		dic['BetaLactamaseAssays_Total'] = row[9]
		dic['FluorescenceAssays_Active'] = row[10]
		dic['FluorescenceAssays_Total'] = row[11]
		return dic
	
	def buildPainsDic(cid, tc):
		cur2 = conn.cursor()
		sqlQuery = "SELECT SMILES, PAINS_A, PAINS_B, PAINS_C, AllAssays_Active, AllAssays_Total, LuciferaseAssays_Active, LuciferaseAssays_Total, BetaLactamaseAssays_Active, BetaLactamaseAssays_Total,"
		sqlQuery += " FluorescenceAssays_Active, FluorescenceAssays_Total, PAINS_A_highlights, PAINS_B_highlights, PAINS_C_highlights FROM painsCompounds where CID = " + str(cid) + ";"
		cur2.execute(sqlQuery)
		row = cur2.fetchone()
		dics = []
		if(row != None):
			#'pains_A' : row[1], 'pains_B' : row[2], 'pains_C': row[3], 'pains_A_highlights' : row[12], 'pains_B_highlights' : row[13], 'pains_C_highlights' : row[14]}
			if((row[1] != None) and (row[12] != None)):	dics.append(rowToDic(cid, tc, row, row[1], row[12]))
			if((row[2] != None) and (row[13] != None)): dics.append(rowToDic(cid, tc, row, row[2], row[13]))
			if((row[3] != None) and (row[14] != None)): dics.append(rowToDic(cid, tc, row, row[3], row[14]))
		return dics

	def buildNonPainsDic(cid,tc):	
		cur2 = conn.cursor()
		sqlQuery = "SELECT SMILES, AllAssays_Active, AllAssays_Total, LuciferaseAssays_Active, LuciferaseAssays_Total, BetaLactamaseAssays_Active, BetaLactamaseAssays_Total, FluorescenceAssays_Active, FluorescenceAssays_Total  FROM nonPainsCompounds WHERE CID= " + str(cid) +" ;"
		#print(sqlQuery)
		cur2.execute(sqlQuery)
		row = cur2.fetchone()
		if(row != None):
			dic = {'cid': cid, 'tc' : tc, 'smile' : row[0]}
			dic['AllAssays_Active'] = row[1]
			dic['AllAssays_Total'] = row[2]
			dic['LuciferaseAssays_Active'] = row[3]
			dic['LuciferaseAssays_Total'] = row[4]
			dic['BetaLactamaseAssays_Active'] = row[5]
                        dic['BetaLactamaseAssays_Total'] = row[6]
                        dic['FluorescenceAssays_Active'] = row[7]
			dic['FluorescenceAssays_Total'] = row[8]
			return dic

		else: return None

	'''
	def findNearest(rows, pains):
		NNs = []
		NN = None
		bestTc = -1.0
		for row in rows:
			fingerprint = row[1]
			compoundFPS = cPickle.loads(str(fingerprint))
			#compoundTc = DataStructs.FingerprintSimilarity(compoundFPS, queryFPS)
			compoundTc = DataStructs.TanimotoSimilarity(compoundFPS, queryFPS)
			if(compoundTc > bestTc): #best match found so far
				if(bestTc > 0.9): NNs.append(NN)
				bestTc = compoundTc
				if(pains):NN = buildPainsDic(row[0],compoundTc)
				else: NN = buildNonPainsDic(row[0],compoundTc)
			elif (compoundTc > 0.9): #not the best match, but still above 0.9
				if(pains):NNs.append(buildPainsDic(row[0],compoundTc))
				else: NNs.append(buildNonPainsDic(row[0],compoundTc))
		return (NN, NNs)
	'''	
	def findNearest(rows, id):
		NNs = []
		for row in rows:
			fingerprint = row[1]
			compoundFPS = cPickle.loads(str(fingerprint))
			compoundTc = DataStructs.TanimotoSimilarity(compoundFPS, queryFPS)
			NNs.append((compoundTc, row[0], id))	
		NN2 = sorted(NNs,cmp=None,key=None, reverse=True)
		return NN2
		
	def getTop5OrAbove90(nearestNeighborTuples, isPains=True): #needs to be sorted before being passed in
		NNs = []
		for i in range(len(nearestNeighborTuples)):
			tup_i = nearestNeighborTuples[i]
			if(i < 5 or tup_i[0] > 0.9):
				if(isPains):
					NNs.extend(buildPainsDic(tup_i[1],tup_i[0]))
				else:
					NNs.append(buildNonPainsDic(tup_i[1],tup_i[0]))
		return NNs
	
#	cur.execute(sqliteQuery1)		
#	rows = cur.fetchall()
#	painsNNsTup = []
#	if(hasPains):
#		painsNNsTup = findNearest(rows, 'painA')
	
	cur.execute(sqliteQuery2)
	rows = cur.fetchall()

	otherPainsNNsTup = findNearest(rows, 'other')

	cur.execute(sqliteQuery3)
	rows = cur.fetchall()

	nonPainsNNsTup = findNearest(rows, 'non')

	painsNNs = []
	otherPainsNNs = []
	nonPainsNNs = []

	if(hasPains):#Pick top 5 for each table. If there are more matches after 5 with hit rate above 0.9, add those too.
		for i in range(len(painsSQLQueries)):
			cur.execute(painsSQLQueries[i])		
			rows = cur.fetchall()
			painsNNsTup = findNearest(rows, 'painA')
			pains_i_NN = getTop5OrAbove90(painsNNsTup)
			painsNNs.append(pains_i_NN)
		otherPainsNNs = getTop5OrAbove90(otherPainsNNsTup)
		nonPainsNNs = getTop5OrAbove90(nonPainsNNsTup, False)
	else:
		nns = []
		neighbor = sorted(otherPainsNNsTup + nonPainsNNsTup)
		for i in range(16):
			nns.append(neighbor.pop())
		nextElem = neighbor.pop()
		while ((nextElem != None) & (nextElem[0] > 0.9)):
			nns.append(nextElem)
			nextElem = neighbor.pop()
		for nn in nns :
			if(nn[2] == 'other'):
				otherPainsNNs.extend(buildPainsDic(nn[1],nn[0]))
			else:
				nonPainsNNs.append(buildNonPainsDic(nn[1],nn[0]))
		
	response = {'pains_nearest' : painsNNs, 'other_pains_nearest' : otherPainsNNs, 'non_pains_nearest':nonPainsNNs, 'flagged_alerts' : flaggedAlerts, 'match_indices': matchIndices}

	
	similar_mols_json = jsonify(response)
	return similar_mols_json
	
	
'''

for row in nonPainsCompoundData:
		compoundMol = Chem.MolFromSmiles(row[1])
		compoundFPS = MACCSkeys.GenMACCSKeys(compoundMol)
		print(compoundFPS)
		compoundTc = DataStructs.FingerprintSimilarity(compoundFPS, queryFPS)
		compoundDic = { 'compound' : row, 'tc' : compoundTc, 'highlights' : highlights }
		if compoundTc > mostSimilarTc : #finding mostSimilarMol
			mostSimilarTc = compoundTc
			if(compoundTc > 0.9): nonPainsSimilarMols.append(nonPainsMostSimilarMol) #append older mostSimilar if above 0.9
		elif compoundTc > 0.9: #Tc value > 0.9
			mostSimilarMol = compoundDic
			nonPainsSimilarMols.append(compoundDic) 
'''	
