from typing import Iterable, Sequence

from .clonemanager import CMLoader
from .sequence import SeqLoaderItem, SeqLoaderManager

def load_sequences(paths: Sequence[str]) -> Iterable[SeqLoaderItem]:
    loader = SeqLoaderManager(
        paths,
        [CMLoader()]
    )

    return loader.load_files()