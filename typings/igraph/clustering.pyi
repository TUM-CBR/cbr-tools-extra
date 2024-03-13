from typing import Iterator, Optional, Sequence

class VertexClustering:

    def __iter__(self) -> Iterator[Sequence[int]]: ...

class VertexDendrogram:
    
    def as_clustering(self, n: Optional[int] = ...) -> VertexClustering: ...