import sqlite3

conn = sqlite3.connect("pains.db")

cur = conn.cursor()

#cur.execute("SELECT * FROM  LIMIT 1")

sqliteQuery2 = "SELECT COUNT(*)  FROM painsCompounds"
cur.execute(sqliteQuery2)
result = cur.fetchone()
print(result,'pains')
sqliteQuery2 = "SELECT COUNT(*)  FROM nonPainsCompounds"
cur.execute(sqliteQuery2)
result = cur.fetchone()
print(result,'non-pains')
'''
	sqliteQuery2 = "SELECT CID, MACCS_fingerprint FROM painsCompounds WHERE NOT (" + wherePains +");"
	
	sqliteQuery3 = "SELECT CID, MACCS_fingerprint FROM nonPainsCompounds;"
'''
