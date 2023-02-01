#!/usr/bin/env python3
import sys
if sys.version_info.major < 3 or sys.version_info.minor < 8:
    sys.exit("gen_ott_properties.py: Python 3.8 or greater required.\n")
import json



major = 3
minor = 4
comp_ext = "tgz"
pub_date = "20221031"
gen_date = "20221031"
amendments_sha = "b6f2e8334447cb2ebe1ff8448cb5872672aaebe9"
gbif_version = "20190916"
genbank_version = "20161216"
idlist_version = "3.1"
irmng_version = "20140131"
ncbi_version = "20191015"
prev_major = "3"
prev_minor = "3"


d = {
  "archive_file_name": f"ott{major}.{minor}.{comp_ext}", 
  "date": f"{pub_date}", 
  "draft": f"1", 
  "generated_on": f"{gen_date}", 
  "legal": f"cc0", 
  "name": f"ott{major}.{minor}", 
  "ott_idspace": f"ott", 
  "separator": f"", 
  "series": f"ott", 
  "sources": {
    "amendments": f"amendments-{amendments_sha}", 
    "fung": f"fung-9", 
    "gbif": f"gbif-{gbif_version}", 
    "genbank": f"genbank-{genbank_version}", 
    "idlist": f"idlist-{idlist_version}", 
    "irmng": f"irmng-{irmng_version}", 
    "ncbi": f"ncbi-{ncbi_version}", 
    "ott": f"ott{prev_major}.{prev_minor}", 
    "silva": f"silva-115", 
    "worms": f"worms-20170604"
  }, 
  "suffix": f".{comp_ext}", 
  "version": f"{major}.{minor}"
}

print(json.dumps(d, sort_keys=True, indent=2))