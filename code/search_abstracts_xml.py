#!/usr/local/bin/python

# Required packages (re, csv, time biopython)
import re
import time
from Bio import Entrez

Entrez.email = 'kipp.johnson@gmail.com'

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

def getxmlAbstract(pmid):
    result = Entrez.efetch(db='pubmed', id=pmid, rettype='abstract', retmode='xml')
    return Entrez.read(result)

def parseAbstract(xmlobj):
    try:
        PMID = str(res[0]['MedlineCitation'].get('PMID').decode())
    except:
        PMID = 'PMID could not be obtained'
    try:
        year = res[0]['MedlineCitation'].get('DateCompleted').get('Year')
    except:
        year = 'Year could not be obtained'
    try:
        author = res[0]['MedlineCitation']['Article']['AuthorList'][0].get('LastName')
    except:
        author = 'Author could not be obtained'
    try:
        journal = res[0]['MedlineCitation']['Article']['Journal']['Title']
    except:
        journal = 'Journal could not be obtained'
    try:
        volume = res[0]['MedlineCitation']['Article']['Journal']['JournalIssue'].get('Volume')
    except:
        volume = 'Volume could not be obtained'
    try:
        issue = res[0]['MedlineCitation']['Article']['Journal']['JournalIssue'].get('Issue')
    except:
        issue = 'Issue could not be obtained'
    try:
        title = res[0]['MedlineCitation']['Article']['ArticleTitle']
    except:
        title = 'Title could not be obtained'
    try:
        abstract = (map(lambda x: ''.join(x) , res[0]['MedlineCitation']['Article']['Abstract']['AbstractText']))
        try:
            map(lambda x: str.encode(x), abstract)
            abstract = " ".join(abstract)
        except:
            print 'Abstract encoding failed, continuing...'
    except:
        abstract = 'Abstract could not be obtained'
    return [PMID, year, author, journal, volume, issue, title, abstract]

if __name__ == "__main__":
    # Get all PMIDs from searchterms
    for term in searchTermList:
        pmidStore.append( getPMIDs(term) )
    # Get only unique PMIDs
    unqPMID = uniquePMIDs(pmidStore)
    print 'Number of unique PMIDs: ', len(unqPMID)

    fo = open("ParsedXMLAbstracts.txt","wb") # in CWD
    fo.write('%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t \n'
    % ("PMID","Year","Author","Journal","Volume","Issue","Title","Abstract"))

    # Fetch abstracts one by one and operate on them
    i = float(0); j=float(len(unqPMID)); # for progress bar

    for pmid in unqPMID:
        res = getxmlAbstract(pmid) # get abstract result from PMID
        ABres = parseAbstract(res) # parse it into usable format
        fo.write('%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t \n' % (ABres[0].encode('utf-8'),
                                                                 ABres[1].encode('utf-8'),
                                                                 ABres[2].encode('utf-8'),
                                                                 ABres[3].encode('utf-8'),
                                                                 ABres[4],
                                                                 ABres[5],
                                                                 ABres[6].encode('utf-8'),
                                                                 ABres[7]))
        time.sleep(0.1) # so the server doesn't reject our requests
        print i/j # progress in console
        i+=1
        if i % 10 == 0: # you can go faster if you sleep the progress occasionally
            time.sleep(5)

        # For testing purposes:
        # if i >19:
        #     break

    # Close the output file
    fo.close()
