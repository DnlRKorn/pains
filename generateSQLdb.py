import sqlite3
import csv
import cPickle
from rdkit import Chem
#from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem

def substruct(compoundMol, painsMol):
	if Chem.AddHs(compoundMol).HasSubstructMatch(painsMol): #search with explicit hydrogens added
		return Chem.AddHs(compoundMol).GetSubstructMatch(painsMol)
		
	elif compoundMol.HasSubstructMatch(painsMol): #search without explicit hydrogens added
		return compoundMol.GetSubstructMatch(painsMol)
		
	return []

conn = sqlite3.connect('pains.db')
#conn.text_factory = str
c = conn.cursor()

#c.execute("DROP TABLE compounds")
c.execute("DROP TABLE painsCompounds")
c.execute("DROP TABLE nonPainsCompounds")

c.execute("CREATE TABLE painsCompounds (CID integer PRIMARY KEY,SMILES varchar(255), PAINS_A varchar(31),PAINS_B varchar(31), PAINS_C varchar(31),AllAssays_Active integer, AllAssays_Total integer,LuciferaseAssays_Active integer, LuciferaseAssays_Total integer, BetaLactamaseAssays_Active integer, BetaLactamaseAssays_Total integer, FluorescenceAssays_Active integer, FluorescenceAssays_Total integer, PAINS_A_highlights string, PAINS_B_highlights string, PAINS_C_highlights string, Mol blob, MACCS_fingerprint blob)")

c.execute("CREATE TABLE nonPainsCompounds (CID integer PRIMARY KEY, SMILES varchar(255), AllAssays_Active integer, AllAssays_Total integer,LuciferaseAssays_Active integer, LuciferaseAssays_Total integer, BetaLactamaseAssays_Active integer, BetaLactamaseAssays_Total integer, FluorescenceAssays_Active integer, FluorescenceAssays_Total integer, Mol blob, MACCS_fingerprint blob)")


conn.commit()

compoundData = []

painsData = {}	#tuples of data about everyPAINS alert taken from master CSV

with open('PAINS_master.csv') as painsAlertCSV:
	reader = csv.DictReader(painsAlertCSV, delimiter=',')
	for row in reader:
		painsData[row['\xef\xbb\xbfALERT NAME']] =  row['SMARTS']
		#painsData.update({row['\xef\xbb\xbfALERT NAME'],row['SMARTS']})

	painsData[' '] = None

counter = 0
	
with open('PAINS_compounds_noDuplicates.csv') as painsCompoundsCSV:
	reader = csv.DictReader(painsCompoundsCSV, delimiter=',')
	for row in reader:
		if(len(row) == 12): continue
#		print row['SMILES']
		compoundMol = Chem.MolFromSmiles(row['SMILES'])
		if(compoundMol != None): 
			painA = row['PAINS_A']
			painB = row['PAINS_B']
			painC = row['PAINS_C']
			if "," in painA:
				print(painA)
				painB = painA.split(',')[1]
				painA = painA.split(',')[0]
			if "," in painB:
				print(painB)
				painA = painB.split(',')[0]
				painB = painB.split(',')[1]
			if "," in painC:
				print(painC)
				painB = painC.split(',')[0]
				painC = painC.split(',')[1]
			molA = Chem.MolFromSmarts(painsData.get(painA))
			molB = Chem.MolFromSmarts(painsData.get(painB))
			molC = Chem.MolFromSmarts(painsData.get(painC))
			painAHighlights = None
			painBHighlights = None
			painCHighlights = None
			if(painA == ' '): painA = None
			if(painB == ' '): painB = None
			if(painC == ' '): painC = None
			
			if(molA != None): 
				hi = substruct(compoundMol,molA)
				painAHighlights = '-'.join([str(i) for i in hi])
			if(molB != None): 
				hi = substruct(compoundMol,molB)
				painBHighlights = '-'.join([str(i) for i in hi])
			if(molC != None): 
				hi = substruct(compoundMol,molC)
				painCHighlights = '-'.join([str(i) for i in hi])
			
			molString = cPickle.dumps(compoundMol)
			macString = cPickle.dumps(AllChem.GetMorganFingerprint(compoundMol, 2))
			
			tempDataTuple =((
			row['\xef\xbb\xbfCID'], 
			row['SMILES'], #2
			painA, 
			painB, 
			painC,#5
			row['AllAssays_Active'], row['AllAssays_Total'], #7 
			row['LuciferaseAssays_Active'], row['LuciferaseAssays_Total'], #9 
			row['BetaLactamaseAssays_Active'], row['BetaLactamaseAssays_Total'], #11 
			row['FluorescenceAssays_Active'], row['FluorescenceAssays_Total'], #13
			painAHighlights,
			painBHighlights,
			painCHighlights, #16
			molString,
			macString # 18
			)) #3 PAINS alerts b/c a compound can have more than 1
			compoundData.append(tempDataTuple)
#			print counter
#			print len(tempDataTuple)
			counter+=1

c.executemany('INSERT INTO painsCompounds VALUES (?,?,?,?,?,?,     ?,?,?,?,?,?,    ?,?,?,?,?,?)', compoundData)
	
conn.commit()

print 'Commited ' + str(len(compoundData)) + ' to table compounds'

compoundData = []

counter = 0

with open('PAINS_Relief_NOPAINS.csv') as painsCompoundsCSV:
	reader = csv.reader(painsCompoundsCSV, delimiter=',')
	for row in reader:
#		print counter
		counter+=1
#		if(len(row) == 9): continue
		#print row['SMILES']
		compoundMol = Chem.MolFromSmiles(row[1])
		if(compoundMol != None):
			molString = cPickle.dumps(compoundMol)
			macString = cPickle.dumps(AllChem.GetMorganFingerprint(compoundMol,2))		
			
			tempDataTuple =((row[0], row[1],
			row[2], row[3], 
			row[4], row[5], 
			row[6], row[7], 
			row[8], row[9],
			molString,
			macString
			))
			compoundData.append(tempDataTuple)
		else: print row[0]

print 'Commited ' + str(len(compoundData)) + ' to table compounds'		
		
c.executemany('INSERT INTO nonPainsCompounds VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', compoundData)
	
print 'Commited ' + str(len(compoundData)) + ' to table compounds'
	
conn.commit()
# We can also close the connection if we are done with it.
# Just be sure any changes have been committed or they will be lost.
conn.close()
