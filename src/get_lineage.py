import urllib, urllib.parse, urllib.request



def get_uniprot_dict(ids, cols):

    # Format the lists of IDs and columns correctly

    cols = ",".join(cols)
    ids = ' or '.join(ids)

    url = 'https://www.uniprot.org/uniprot/'

    params = {
        'format': 'tab',
        'query': ids,
        'columns': "id," + cols
    }

    data = urllib.parse.urlencode(params).encode('utf-8')
    request = urllib.request.Request(url, data)
    print (data)
    opener = urllib.request.build_opener()
    response = opener.open(request)
    page = response.read(200000).decode('utf-8')

    up_dict = {}

    # For each record we retrieve, split the line by tabs and build up the UniProt dict
    for line in page.split("\n")[1:]:
        if line:
            splitlines= line.split("\t")
            id_dict = {}
            pos = 1
            for col in cols.split(","):
                id_dict[col] =splitlines[pos]
                pos +=1
            up_dict[splitlines[0]] = id_dict

    return up_dict


# Usage - with a list of UniProt IDs

up_dict = get_uniprot_dict(["P05791", "Q9LIR4"], ["lineage(KINGDOM)","lineage(SPECIES)"])


print (up_dict)
