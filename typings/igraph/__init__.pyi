from igraph import clustering
from igraph.seq import VertexSeq
from pandas import DataFrame
from typing import Any, Literal, Optional, Sequence, Union

GraphSearchMode = Union[Literal['strong'], Literal['weak']]

class Vertex:

    def __getitem__(self, key: str) -> Any: ...


EdgeWeights = Sequence[float]

class Graph:

    def connected_components(self,  mode: GraphSearchMode = ...) -> clustering.VertexClustering: ...

    @property
    def vs(self) -> VertexSeq: ...

    def community_fastgreedy(self, weights: Optional[EdgeWeights] = ...) -> clustering.VertexDendrogram: ...

    def community_walktrap(self, weights: Optional[EdgeWeights] = ..., steps: int = ...) -> clustering.VertexDendrogram: ...

    @classmethod
    def DataFrame(
        cls,
        edges: DataFrame,
        directed: bool = ...,
        vertices: Optional[DataFrame] = ...,
        use_vids: bool = ...,
    ) -> 'Graph': ...