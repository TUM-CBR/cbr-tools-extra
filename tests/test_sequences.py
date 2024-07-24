import pytest

from cbrextra.sequences.data import InteractiveInput, InteractiveOutput, SearchArg, SearchArgs
from cbrextra.sequences.search import Search
from cbrextra.sequences.main import SequencesInteractive
from cbrextra import testutils
from os import path

from .support import *

class TestSequences:

    def test_search(self):

        expected = {
            "test_search": [
                "GCA_008369605.1",
                "GCA_000969265.1",

            ],
            "test_search_2": [
                "GCA_002634065.1"
            ]
        }
        search = Search()
        search_args = SearchArgs(
            searches=[
                SearchArg(
                    search_id="test_search",
                    tax_ids=[666],
                ),
                SearchArg(
                    search_id="test_search_2",
                    names=["Coxiella burnetii"],
                )
            ]
        )

        results = list(search.run(search_args))

        for result in results:
            search_id = result.search_id
            accessions = { r.accession for r in result.records }

            for exp_value in expected[search_id]:
                assert exp_value in accessions, f"Accessio {expected} was not present" 

    def test_interactive(
        self,
        tmp_folder
    ):

        db_file = path.join(tmp_folder, "sequences.db")
        search = SequencesInteractive(db_file=db_file)

        in_messages = [
            InteractiveInput(
                search=SearchArgs(
                    searches = [
                        SearchArg(
                            search_id="1",
                            tax_ids=[666]
                        ),
                        SearchArg(
                            search_id="2",
                            names=["Coxiella burnetii"]
                        )
                    ]
                )
            )
        ]

        for result in testutils.test_interactive(search, in_messages, InteractiveInput, InteractiveOutput):
            print(result)
            value = result.result
            assert value is not None
            assert value.model_dump_json() is not None
