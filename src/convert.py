#!/usr/bin/python

import dendropy
from datetime import datetime
import re

pedv = dendropy.DnaCharacterMatrix.get(path="pedv.fasta", schema="fasta")

for taxa in pedv.taxon_namespace:
    date = re.search("([0-9]{2}\-[A-Z][a-z]{2}\-[0-9]{4})", taxa.label).group(1)
    date = datetime.strptime(date, '%d-%b-%Y').date()
    decdate = date.year + (date.month - 1) / 12  + date.day / 365
    decdate = "{0:.2f}".format(decdate)
    taxa.label = taxa.label + '_' + decdate

pedv.write(path="tipdate.in", schema="phylip")
