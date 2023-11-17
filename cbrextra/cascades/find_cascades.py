from .store import Store

def find_cascades(
    results_db : str,
    stop_identity_treshold : float
):
    store = Store(results_db)