from igraph.seq import VertexSeq
from pandas import DataFrame
from typing import Any, Iterator, List, Literal, Optional, Union

class VertexClustering:

    def __iter__(self) -> Iterator[List[int]]: ...

GraphSearchMode = Union[Literal['strong'], Literal['weak']]

class Vertex:

    def __getitem__(self, key: str) -> Any: ...

class Graph:

    def connected_components(self,  mode: GraphSearchMode = ...) -> VertexClustering: ...

    @property
    def vs(self) -> VertexSeq: ...

    @classmethod
    def DataFrame(
        cls,
        edges: DataFrame,
        directed: bool = ...,
        vertices: Optional[DataFrame] = ...,
        use_vids: bool = ...,
    ) -> 'Graph': ...