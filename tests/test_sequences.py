from cbrextra.sequences.data import SearchArg, SearchArgs
from cbrextra.sequences.search import Search
from cbrextra import testutils

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
        
