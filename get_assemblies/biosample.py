#!/usr/bin/env python
from rich.console import Console
from rich.pretty import pprint

console = Console(stderr=True)
KEYS = [
"env_broad_scale",
"env_local_scale",
"geo_loc_name",
"env_medium",
"depth",
"isolation source",
"lat_lon",
"collection_date",
]

class BioSample:
    def __init__(self, acc: str, assemid: str, attrs: list = KEYS):
        self.acc = acc
        self.assemid = assemid
        self.attrs = attrs
        self.isolate = '-'
        self.initattrs(attrs)

    def initattrs(self, attrs:list):
        for attr in attrs:
            setattr(self, attr, '-')

    def __rich_repr__(self):
        yield "acc", self.acc
        yield "assemid", self.assemid
        yield "isolate", self.isolate
        for attr in self.attrs:
            yield attr, getattr(self, attr)

    def pprint(self):
        pprint(self, console=console)

def get_bs2assemid(dl_mapping: dict):
    """
    dl_mapping keys are assembly accessions.
    "biosample" is a key of the dl_mapping[acc] dict.
    """
    bs2assemid = {}
    for acc in dl_mapping:
        bsid = dl_mapping[acc]["biosample"]
        bs2assemid[bsid] = acc
    return bs2assemid
