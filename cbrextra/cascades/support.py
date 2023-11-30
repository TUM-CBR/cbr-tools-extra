from Bio.Blast.NCBIWWW import QBlastStatusMonitor
from concurrent.futures import Executor
import json
from typing import List, TextIO

from .store import FindOrganismsArgs, Iterable, NamedTuple, Optional, Store

class RunCascadesContext(NamedTuple):
    store : Store
    executor : Executor
    out_stream : TextIO
    domain : Optional[str]

    def write_status(self, info: dict):
        info['type'] = 'status'

        self.out_stream.write(json.dumps(info))
        self.out_stream.write("\n")
        self.out_stream.flush()

    def write_error(
        self,
        step: int,
        error : Exception,
        extra : Optional[dict] = None
    ):
        extra = extra or {}
        extra['step'] = step
        extra['error'] = str(error)

        self.write_status(extra)

    def with_domain(self, arg : FindOrganismsArgs) -> FindOrganismsArgs:
        if self.domain is None:
            return arg

        return arg._replace(included_organisms=[self.domain])

    def with_domains(self, args: Iterable[FindOrganismsArgs]) -> List[FindOrganismsArgs]:
        return [self.with_domain(arg) for arg in args]


class BlastMonitor(QBlastStatusMonitor):

    def __init__(
        self,
        step : int,
        context : RunCascadesContext
    ):
        self.__step = step
        self.__context = context

    def on_status(self, payload: str) -> None:
        self.__context.write_status({
            'message': payload,
            'step': self.__step
        })

    def on_timeout(self, e: TimeoutError) -> None:
        self.__context.write_error(
            self.__step,
            e
        )

    def set_status_request(self, url: str, msg: bytes) -> None:
        self.__context.write_status({
            'url': url,
            'msg': msg.decode(),
            'step': self.__step
        })