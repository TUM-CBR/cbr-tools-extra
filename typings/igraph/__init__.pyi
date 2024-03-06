from igraph.seq import VertexSeq
from typing import Iterator, List, Literal, Union

class VertexClustering:

    def __iter__(self) -> Iterator[List[int]]: ...

GraphSearchMode = Union[Literal['strong'], Literal['weak']]

class Graph:

    def connected_components(self,  mode: GraphSearchMode = ...) -> VertexClustering: ...

    @property
    def vs(self) -> VertexSeq: ...