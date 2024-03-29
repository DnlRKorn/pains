from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys, AllChem, rdmolops
import csv
from flask import jsonify
import sqlite3
import cPickle

########## Populating arrays from CSV files ########## 

def main():
	smileStr=raw_input()
	pains(smileStr)
	
	
def pains(smileStr):
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
		#	for i in range(len(all_matches)):
		#		b = True
		#		for j in range(i):
		#			diff = 0
		#			for x, y in zip(all_matches[i], all_matches[j]):
		#				diff = diff + abs(x - y)
			#		if(diff > 12): b = False #diff is very large, two highlights are similar
			#	if(b): matches.append(all_matches[i])
			matches.append(all_matches[0])
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
		tmp_str = "PAINS_A = '" + str(match[1]) + "' OR PAINS_B = '"+ str(match[1]) + "' OR PAINS_C = '"+ str(match[1]) + "' "
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
	
	
	sqliteQuery2 = "SELECT CID, MACCS_fingerprint FROM painsCompounds WHERE NOT (" + wherePains +");"
	
	sqliteQuery3 = "SELECT CID, MACCS_fingerprint FROM nonPainsCompounds;"
	
	conn = sqlite3.connect("pains.db")
	#conn.text_factory = str
	
	cur = conn.cursor()
	cur2 = conn.cursor()
	queryFPS = AllChem.GetMorganFingerprint(queryMol,2)

	
	def rowToDic(cid, tc, row, pains, pains_highlights):
		dic = {'cid': cid, 'tc' : tc, 'smile' : row[0]}
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
		print('Pains A ' + str(row[1]) + ' , B ' + str(row[2]) + ' , C ' + str(row[3]) + ' highlight A ' + str(row[12]) + ' , highlight B ' + str(row[13]) + ' , highlight C ' + str(row[14]))
		if(row != None):
			#'pains_A' : row[1], 'pains_B' : row[2], 'pains_C': row[3], 'pains_A_highlights' : row[12], 'pains_B_highlights' : row[13], 'pains_C_highlights' : row[14]}
			if((row[1] != None) and (row[12] != None) and ('-' in row[12])): dics.append(rowToDic(cid, tc, row, row[1], row[12]))
			if((row[2] != None) and (row[13] != None) and ('-' in row[13])): dics.append(rowToDic(cid, tc, row, row[2], row[13]))
			if((row[3] != None) and (row[14] != None) and ('-' in row[14])): dics.append(rowToDic(cid, tc, row, row[3], row[14]))
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

	#Takes a list of compounds and generates a list of tuples sorted by Tc similarity to the query mol.
	#Allows a unique identifier to be passed in and appended at the end of each tuple.
	def findNearest(cur, sqliteQuery, id):
		tup = (-1, -1, id) #essentially negitive similarity
		NNs = []
		for i in range(20): NNs.append(tup)
		for row in cur.execute(sqliteQuery):
			fingerprint = row[1]
			compoundFPS = cPickle.loads(str(fingerprint))
			compoundTc = DataStructs.TanimotoSimilarity(compoundFPS, queryFPS)
			for i in range(20):
				if(compoundTc > NNs[i][0]):
					NNs.insert(i, (compoundTc, row[0], id))
					NNs.pop()
					break
		#sortedNN = sorted(NNs,cmp=None,key=None, reverse=True)
		print('len of sorted ' + str(len(NNs)) + " " + str(id))
		return NNs
	
	#Grabs all enteries in the list of tuples. Return all those that are above 90% Tc Similarity.
        #Return top 5 no matter what their similarity is.	
	def getTop5OrAbove90(nearestNeighborTuples, isPains=True): #needs to be sorted before being passed in
		NNs = []
		#print('len of nn tuples ' + str(len(nearestNeighborTuples)))
		if(isPains):
			for i in range(len(nearestNeighborTuples)):
				tup_i = nearestNeighborTuples[i]
				print(tup_i)
				if(len(NNs) < 5 or tup_i[0] > 0.9):
					NNs.extend(buildPainsDic(tup_i[1],tup_i[0]))
		else:
			for i in range(len(nearestNeighborTuples)):
				tup_i = nearestNeighborTuples[i]
				if(len(NNs) < 5 or tup_i[0] > 0.9):
					NNs.append(buildNonPainsDic(tup_i[1],tup_i[0]))

		return NNs
	
	#cur.execute(sqliteQuery2)
	#rows = cur.fetchall()
	#print(len(rows))
	#print('others lenght')
	otherPainsNNsTup = findNearest(cur, sqliteQuery2, 2)

	cur.execute(sqliteQuery3)
	rows = cur.fetchall()

	nonPainsNNsTup = findNearest(cur, sqliteQuery3, 0)

	painsNNs = []
	otherPainsNNs = []
	nonPainsNNs = []
	print('has pains ' + str(hasPains))
	if(hasPains):#Pick top 5 for each table. If there are more matches after 5 with hit rate above 0.9, add those too.
		for i in range(len(painsSQLQueries)):
			print(painsSQLQueries[i])
			#cur.execute(painsSQLQueries[i])		
			#rows = cur.fetchall()
			painsNNsTup = findNearest(cur, painsSQLQueries[i], 'painA')
			pains_i_NN = getTop5OrAbove90(painsNNsTup)
			painsNNs.append(pains_i_NN)
		print('other pains tup ' + str(len(otherPainsNNsTup)))
		otherPainsNNs = getTop5OrAbove90(otherPainsNNsTup)
		print('other pains nn ' + str(len(otherPainsNNs)))
		nonPainsNNs = getTop5OrAbove90(nonPainsNNsTup, False)

		print('Non pains nn ' + str(len(nonPainsNNs)))

	else:
		nns = []
		neighbor = sorted(otherPainsNNsTup + nonPainsNNsTup)
		for i in range(16):
			nns.append(neighbor.pop())
		nextElem = neighbor.pop()
		while ((nextElem != None) & (nextElem[0] > 0.9)):
			nns.append(nextElem)
			nextElem = neighbor.pop()
		print(nns[0])
		nextElem = neighbor.pop()
		while (len(otherPainsNNs) + len(nonPainsNNs)) < 15 or nextElem[0] > 0.9:
			nn = nextElem
			if(nn[2] == 2):
				otherPainsNNs.extend(buildPainsDic(nn[1],nn[0]))
			else:
				nonPainsNNs.append(buildNonPainsDic(nn[1],nn[0]))
			nextElem = neighbor.pop()
		
	response = {'pains_nearest' : painsNNs, 'other_pains_nearest' : otherPainsNNs, 'non_pains_nearest':nonPainsNNs, 'flagged_alerts' : flaggedAlerts, 'match_indices': matchIndices}

	
	similar_mols_json = jsonify(response)
	return similar_mols_json
