"""Usage:
    cbrtools sequences scan --db-file=<db-file> [-- <scan-dir>...]
"""

from typing import Optional, Sequence
from cbrextra.core.module import Context, Result
from docopt import docopt

from ..core.module import Module
from .procedures import run_scan

class SequencesModule(Module):

    def scan(self, db_file: str, scan_dirs: Optional[Sequence[str]]):

        scan_dirs = [] if scan_dirs is None else scan_dirs

        run_scan(db_file, scan_dirs)

        return Result.success()

    def main(self, context: Context) -> Result:

        options = docopt(__doc__)

        if options.get("scan"):
            return self.scan(options['--db-file'], options.get('scan-dir'))
        
        return Result.not_requested()