import urllib, urllib.parse, urllib.request



def map_ids():

    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from':'PDB_ID',
    'to':'ID',
    'format':'tab',
    'query':'3QM4 '
    }

    data = urllib.parse.urlencode(params).encode('utf-8')
    print (data)
    request = urllib.request.Request(url, data)
    print (request)
    contact = "gabriel.foley@uqconnect.edu.au"
    request.add_header('User-Agent', 'Python %s' % contact)

    opener = urllib.request.build_opener()
    response = opener.open(request)

    page = response.read(200000).decode('utf-8')



    print (page)

def query_uniprot():
    pass


map_ids()
