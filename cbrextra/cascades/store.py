from os import path
from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from .models import Base

class Store:

    class StoreApi:

        def __init__(self, session: Session):
            self.__session = session

        def __enter__(self):
            pass

        def __exit__(self, *args, **kwargs):
            self.__session.__exit__(*args, **kwargs)

    def __init__(self, db_file: str):
        self.__db_file = db_file
        self.__engine = create_engine(f"sqlite:///{self.__db_file}")
        self.__init_db()

    def __init_db(self):
        if path.exists(self.__db_file):
            return

        Base.metadata.create_all(self.__engine)

    def __session_factory(self) -> sessionmaker[Session]:
        return sessionmaker(bind = self.__engine)

    