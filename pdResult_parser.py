import sqlite3

def fetch_sqlite3_db(infile, table):

    con = sqlite3.connect(infile)
    cur = con.cursor()
    if table == 'target':
        cur.execute("SELECT * FROM TargetPsms")
        names = list(map(lambda x: x[0], cur.description))
        return names, cur.fetchall()
        
    elif table == 'decoy':
        cur.execute("SELECT * FROM DecoyPsms")
        names = list(map(lambda x: x[0], cur.description))
        return names, cur.fetchall()

    con.close()

