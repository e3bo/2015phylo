#!/usr/bin/python2

from Bio import SeqIO

gb=SeqIO.parse('pedv-north-america.gb','genbank')

recs = []
for g in gb:
    for feat in g.features:
        if feat.type == 'source':
            if 'collection_date' in feat.qualifiers and 'country' in feat.qualifiers:
                rec =  [str(g.id), feat.qualifiers['country'][0]]
                rec.append(feat.qualifiers['collection_date'][0])
                if 'strain' in feat.qualifiers:
                    rec.append(feat.qualifiers['strain'][0])
                rec = ','.join(rec)
                recs.append(rec)
            
f = open('pedv-north-america.csv', 'w')
f.write('country,date,strain\n')
for rec in recs:
    f.write(rec + '\n')
f.close()
