import igraph
from typing import Iterator, List

VertexSeqIndex = List[int]

class VertexSeq:

    def __getitem__(self, ix: VertexSeqIndex) -> 'VertexSeq': ...

    def __iter__(self) -> Iterator['igraph.Vertex']: ...

    def subgraph(self) -> 'igraph.Graph': ...