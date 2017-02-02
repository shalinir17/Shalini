import requests, sys
 
server = "http://grch37.rest.ensembl.org"
ext = "/lookup/id/ENSG00000157764?expand=1"
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()
print repr(decoded)