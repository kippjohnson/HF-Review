from Bio import Entrez

Entrez.email = 'kipp.johnson@icahn.mssm.edu'

search1 = 'congestive heart failure[MeSH Terms]'
result1 = Entrez.esearch(db='pubmed', retmax=100000000, term=search1)
record1 = Entrez.read(result1)
pmids1 = list(record1['IdList'])

search2 = '(patient readmission[MeSH Terms])'#' AND readmi$ rehosp$'
result2 = Entrez.esearch(db='pubmed', retmax=100000000, term=search2)
record2 = Entrez.read(result2)
pmids2 = list(record2['IdList'])

search3 = '(risk[MeSH Terms])'#' AND model$ predict$ use$ util$ risk$'
result3 = Entrez.esearch(db='pubmed', retmax=100000000, term=search3)
record3 = Entrez.read(result3)
pmids3 = list(record3['IdList'])

# print "PMID search 1: ", len(pmids1)
# print "PMID search 2: ", len(pmids2)
# print "PMID search 3: ", len(pmids3)

pmids1 = set(pmids1)
pmids2 = set(pmids2)
pmids3 = set(pmids3)

pmidSet = set.intersection(pmids1, pmids2, pmids3)

print "PMID intersection: ", len(pmidSet)

# Convert from set to list data structure
pmidSet = list(pmidSet)

for pmid in pmidSet:
    res = Entrez.efetch(db='pubmed', id=pmid, rettype='abstract', retmode='text')
    print res.read(), "\n"
