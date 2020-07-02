from suds.client import Client
import suds
import time

import fitlib.helplib as hl

#be aware that you need a chemspider_token.txt in the directory for the app to work
#the chemspider_token.txt should only contain the token (available online for free)

class ChemicalObject():
    def __init__(self, name = '', cas = '', inchi = '', inchikey = '', csid = ''):

        #first define the SOAP service for searching
        searchurl = 'http://www.chemspider.com/Search.asmx?WSDL'

        try:
            self.searchclient = Client(searchurl)
        except Exception as e:
            print(e)

        #define the soap service for inchi-conversion
        inchiurl = 'http://www.chemspider.com/InChI.asmx?WSDL'

        try:
            self.inchiclient = Client(inchiurl)
        except Exception as e:
            print(e)

        #set all properties to the ones from initiating call
        self.cas = cas
        self.inchi = inchi
        self.inchikey = inchikey
        self.name = name
        self.csid = csid

        #no transaction id for now
        self.transaction_id = ''
        #how quickly should we ask for results? in seconds
        self.timetick = 0.2
        self.maxtime = 15
        
        #read chemspider token from config file 'chemspider_token'
        try:
            f = hl.openfile('chemspider_token.txt')
        except IOError:
            raise IOError

        self.token = f.readline()

    def complete(self):
        """Fills all other properties of an instance"""
        #first retrieve the Chemspider ID
        self.retrieve_csid()#
        #fill up the other fields
        self.fill_forms_with_csid()

    def status(self):
        #if we don't have a transaction id, we are free
        if self.transaction_id != '':
            return 'busy'
        else:
            return 'free'

    def retrieve_csid(self):
        #for what should we search?
        if self.inchi != '':
            searchterm = self.inchi
        else:
            searchterm = self.name
            
        #it's a good idea to only search for ascii:
        searchterm = searchterm.decode('utf8', 'replace').encode('ascii', 'replace')
        
        #try connecting
        try:
            self.transaction_id = self.searchclient.service.AsyncSimpleSearch(searchterm, self.token)
        except suds.WebFault as detail:
            self.errorstatus = detail
            self.transaction_id = ''
        
        #don't run too long in the following loop
        i = 0

        #if successful we can check whether the results are already available
        if self.transaction_id != '':
            while self.searchclient.service.GetAsyncSearchStatus(self.transaction_id, self.token) != 'ResultReady':
                #wait a little
                time.sleep(self.timetick)
                i = i + 1
                if i > (self.maxtime / self.timetick):
                    print('No result, aborting')
                    break
                
            #ready! the [0] is because it basically gives a list and we use the first one
            result = self.searchclient.service.GetAsyncSearchResult(self.transaction_id, self.token)
            if result != '':
                self.csid = result[0][0]
        #transaction over. set transaction id to empty for proper status displays
        self.transaction_id = ''

    def fill_forms_with_csid(self):
        """Retrieve all data from Chemspider service using a CS ID"""
        if self.csid != '':
            try:
                tmp = self.searchclient.service.GetCompoundInfo(self.csid, self.token)
            except suds.WebFault as detail:
                print(detail)
            self.inchi = tmp[1]
            self.inchikey = tmp[2]
