#!/usr/bin/env python3
"""
synapse_io.py -- minimal synapseclient helpers for the Python plotting scripts.

Python equivalent of the lab's R synExtra workflow:

    synapser::synLogin()                        ->  login()
    syn("synXXXX")                              ->  fetch("synXXXX")         (download -> local path)
    synStoreMany(paths, parentId=, used=)       ->  store(path, parent, used=, executed=)

AUTH: synapseclient reads credentials from the SYNAPSE_AUTH_TOKEN environment
variable or ~/.synapseConfig -- exactly like `synapse login` / synapser::synLogin().
No token is ever written into the code.
"""
import functools

import synapseclient
from synapseclient import File


@functools.lru_cache(maxsize=1)
def login():
    """Return a cached, logged-in Synapse client (personal-access-token auth)."""
    syn = synapseclient.Synapse()
    syn.login()                 # SYNAPSE_AUTH_TOKEN env var, or ~/.synapseConfig
    return syn


def fetch(syn_id, download_location=None):
    """Download an entity by synID and return its local path (like synExtra's syn())."""
    entity = login().get(syn_id, downloadLocation=download_location)
    return entity.path


def store(path, parent_id, used=None, executed=None, force_version=False):
    """Upload `path` into a Synapse folder with provenance (like synStoreMany).

    used     : list of input synIDs / URLs this output was derived from
    executed : the script that produced it (path or GitHub URL)
    """
    entity = File(path, parent=parent_id)
    return login().store(entity, used=used, executed=executed,
                         forceVersion=force_version)
