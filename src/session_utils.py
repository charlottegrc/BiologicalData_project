from __future__ import annotations
import re
from typing import Dict, Optional
import requests
from requests.adapters import HTTPAdapter, Retry

UNIPROT_REST = "https://rest.uniprot.org"
EBI_PROTEINS_API = "https://www.ebi.ac.uk/proteins/api"
EBI_TAXONOMY = f"{EBI_PROTEINS_API}/taxonomy"

NS = {"uniprot": "http://uniprot.org/uniprot"}  # from api_and_parsing.py pattern

def make_session() -> requests.Session:
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session

# Handle UniProt REST pagination.
def get_next_link(headers: Dict[str, str]) -> Optional[str]:
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)
    return None


def get_batch(session: requests.Session, batch_url: str):
    """Yield (response, total_results) for UniProt REST paginated endpoints."""
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers.get("x-total-results")
        yield response, total
        batch_url = get_next_link(response.headers)

re_next_link = re.compile(r'<(.+)>; rel="next"')