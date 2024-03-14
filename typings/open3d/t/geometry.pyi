from typing import Any
import open3d.core as core
import open3d.geometry as geometry

class RaycastingScene:
    def add_triangles(self, mesh: 'TriangleMesh') -> Any: ...

    def compute_signed_distance(self, query: core.Tensor) -> core.Tensor: ...

class TriangleMesh:

    @classmethod
    def from_legacy(cls, mesh: geometry.TriangleMesh) -> 'TriangleMesh': ...