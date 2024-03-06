import igraph
from typing import List

VertexSeqIndex = List[int]

class VertexSeq:

    def __getitem__(self, ix: VertexSeqIndex) -> 'VertexSeq': ...

    def subgraph(self) -> 'igraph.Graph': ...