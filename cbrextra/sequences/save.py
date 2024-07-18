from typing import Iterable, List

from .data import SearchResultRecord
from .manager import SessionManager

class Save:

    def __init__(
        self,
        db_file: str
    ):
        self.__session_manager = SessionManager.from_file(db_file)

    def save_search_results(
        self,
        records: Iterable[SearchResultRecord]
    ) -> List[Exception]:

        with self.__session_manager.session() as session:
            return session.save_search_results(records)
