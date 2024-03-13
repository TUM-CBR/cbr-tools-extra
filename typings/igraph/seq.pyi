import igraph
from typing import Iterator, Sequence

VertexSeqIndex = Sequence[int]

class VertexSeq():

    def __len__(self) -> int: ...

    def __getitem__(self, ix: VertexSeqIndex) -> 'VertexSeq': ...

    def __iter__(self) -> Iterator['igraph.Vertex']: ...

    def subgraph(self) -> 'igraph.Graph': ...