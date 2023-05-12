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

def fetch_massSpec_info(infile):
    con = sqlite3.connect(infile)
    cur = con.cursor()
    cur.execute("SELECT * FROM MassSpectrumItems")
    names = list(map(lambda x: x[0], cur.description))

    photo = cur.fetchone()[-1]
    print (photo)
    for i in cur.fetchall():
        print (i)
