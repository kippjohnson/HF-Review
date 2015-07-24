#!/usr/local/bin/python

# Required packages (re, csv, time biopython)
import re
import csv
import time
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
    return result.read()

def parseAbstract(abstract):
    result = filter(None, map(lambda i: re.sub("\n"," ", i), abstract.split("\n\n")))
    pmid = result[-1].lstrip("PMID: ").rstrip("  [PubMed - indexed for MEDLINE]")
    journal = result[0].lstrip('1. ')
    title = result[1]
    authors = result[2]
    #author_info = result[3]
    #abstracttext = result[4]
    return([pmid, journal, title, authors])

if __name__ == "__main__":
    # Get all PMIDs from searchterms
    for term in searchTermList:
        pmidStore.append( getPMIDs(term) )
    # Get only unique PMIDs
    unqPMID = uniquePMIDs(pmidStore)
    print 'Number of unique PMIDs: ', len(unqPMID)

    fo = open("ParsedAbstracts.csv","wb") # in CWD
    fo.write('%s\t%s\t%s\t%s\t \n' % ("PMID","Journal","Title","Authors"))

    # Fetch abstracts one by one and operate on them
    i = float(0); j=float(len(unqPMID))
    for pmid in unqPMID:
        res = getAbstract(pmid) # get abstract result from PMID
        ABres = parseAbstract(res) # parse it into usable format
        fo.write('%s\t%s\t%s\t%s\t \n' % (str(ABres[0]),str(ABres[1]),str(ABres[2]),str(ABres[3])) )
        time.sleep(3) # so the server doesn't reject our requests
        print i/j # progress in console
        i+=1

    # Close the output file
    fo.close()
