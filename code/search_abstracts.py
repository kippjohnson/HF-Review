from Bio import Entrez
Entrez.email = 'kipp.johnson@icahn.mssm.edu'

searchTermList = ['congestive heart failure[MeSH Terms]',
                  'patient readmission[MeSH Terms]',
                  'risk[MeSH Terms]']

pmidStore = list()

def getPMIDs(searchterm):
    result = Entrez.esearch(db='pubmed', retmax=10000000, term=searchterm)
    record = Entrez.read(result)
    pmids = list(record['IdList'])
    return pmids

def uniquePMIDs(pmidlist):
    pmid = reduce(set.intersection, map(set, pmidlist))
    pmid = list(pmid)
    return pmid

def getAbstract(pmid):
    result = Entrez.efetch(db='pubmed', id=pmid, rettype='abstract', retmode='text')
    return result

if __name__ == "__main__":
    # Get all PMIDs from searchterms
    for term in searchTermList:
        pmidStore.append( getPMIDs(term) )
    # Get only unique PMIDs
    unqPMID = uniquePMIDs(pmidStore)
    # Fetch abstracts 1 by one and print them
    for pmid in unqPMID:
        res = getAbstract(pmid)
        print(res.read())
