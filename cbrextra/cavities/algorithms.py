from igraph import Graph
import numpy as np
from numpy import float64
from numpy.typing import NDArray
from open3d.geometry import PointCloud, Octree, OctreeInternalPointNode, OctreeNode, OctreeNodeInfo
from open3d.utility import Vector3dVector
import pandas as pd
from typing import cast, Dict, List, NamedTuple

from .data import *

class FindCavitiesResult(NamedTuple):
    context: 'FindCavitiesContext'
    graphs: List[FindCavitiesGraph]
    nodes_df: pd.DataFrame
    edges_df: pd.DataFrame
    depths_to_volume: Dict[int, float]

    def __get_boxes(self, graph: Graph) -> ProteinCavity:
        node_ids = [node['name'] for node in graph.vs]
        nodes_df = self.nodes_df.loc[node_ids]

        return ProteinCavity(nodes_df)

    def __get_cavities_for(self, cavities: FindCavitiesGraph, min_volume: int, max_volume: int):

        graph = cavities.graph
        unit_volume = self.depths_to_volume[cavities.depth]
        min_units = min_volume / unit_volume
        max_units = max_volume / unit_volume

        groups = graph.connected_components()
        accepted = [
            self.__get_boxes(graph.vs[group].subgraph())
            for group in groups
            if len(group) <= max_units and len(group) >= min_units
        ]

        return accepted

    def get_cavities(self, min_volume: int = 2, max_volume: int = 1000) -> List[ProteinCavity]:
        values = [
            cavity
            for graph in self.graphs
            for cavity in self.__get_cavities_for(graph, min_volume, max_volume)
        ]
        return values

class FindCavitiesContext(NamedTuple):
    points: NDArray[float64]
    point_cloud: PointCloud
    octree: Octree
    empty_region_treshold: int

    @staticmethod
    def construct(points: NDArray[float64], empty_region_treshold: int) -> 'FindCavitiesContext':

        pcd = PointCloud()
        pcd.points = Vector3dVector(points)

        # For protein structures, this leads to the smallest boxes
        # measuring almost 1A. This is sufficient for most cases.
        octree = Octree(max_depth=6)
        octree.convert_from_point_cloud(pcd)

        return FindCavitiesContext(
            points=points,
            point_cloud=pcd,
            octree=octree,
            empty_region_treshold=empty_region_treshold
        )

    def get_empty_corners(self) -> pd.DataFrame:
        """Constructs a DataFrame that contains the coordinates of the corners of all of the
        nodes of the Octree which are empty."""

        corners: List[NDArray[float64]] = []
        N_CORNERS = 8

        # Enumerate children coordinates just like
        # Open3d does the enumeration
        # https://github.com/isl-org/Open3D/blob/f5f672b4af1fc81e423c3c1b6215497f5a8816ea/cpp/open3d/geometry/Octree.cpp#L700
        children_coords_offset = np.array([
            [i % 2, int(i/2) % 2, int(i/4) % 2]
            for i in range(N_CORNERS)
        ])

        to_center = np.array([1,1,1])
        box_id = 0

        def apply(node: OctreeNode, node_info: OctreeNodeInfo):
            nonlocal corners, box_id

            child_size = node_info.size / 2

            if isinstance(node, OctreeInternalPointNode):
                for i,child in enumerate(node.children):
                    if child is not None and len(child.indices) > self.empty_region_treshold:
                        continue

                    box_id += 1
                    uid = box_id
                    child_origin = node_info.origin + children_coords_offset[i]*child_size
                    center = child_origin + child_size/2 * to_center

                    # Convert corner's floating point coordinates as an int. We only need the first
                    # decimal position as that gives you a resolution of 0.1A, more than enough for our
                    # purposes. Representing the positions as integers allows one to quickly align corners
                    # using data manipulation libraries
                    box_vertices = ((np.tile(child_origin, (N_CORNERS, 1)) + children_coords_offset*child_size)*10).astype('int')

                    # We create an array with the fields common to every record, then we repeated
                    # once for every vertex of the box (8 times in 3d space)
                    common_values = np.repeat(
                        [
                            np.concatenate([
                                [uid, child_size, node_info.depth + 1],
                                center
                            ])
                        ],
                        len(box_vertices), axis=0
                    )

                    corners.append(np.hstack([common_values, box_vertices]))


        self.octree.traverse(apply)
        values = np.concatenate(corners)

        return pd.DataFrame({
            K_BOX_ID: values[:,0].astype('int'),
            K_BOX_SIZE: values[:,1],
            K_BOX_DEPTH: values[:,2].astype('int'),
            K_BOX_X: values[:,3],
            K_BOX_Y: values[:,4],
            K_BOX_Z: values[:,5],
            K_BOX_CX: values[:,6].astype('int'),
            K_BOX_CY: values[:,7].astype('int'),
            K_BOX_CZ: values[:,8].astype('int')
        })
    
    def construct_graph(self, corners: pd.DataFrame) -> FindCavitiesResult:

        CORNER_JOIN_KEYS = [K_BOX_CX, K_BOX_CY, K_BOX_CZ, K_BOX_DEPTH]
        K_BOX_ID_FROM = "box_id_from"
        K_BOX_ID_TO = "box_id_to"

        # Match the corners of each box with corners of other boxes
        # of the same depth. As all boxes come from an octree, adjecent
        # boxes will always share corners
        aligned = corners.rename({K_BOX_ID: K_BOX_ID_FROM}, axis=1) \
            .merge(
                corners.rename({K_BOX_ID: K_BOX_ID_TO}, axis=1),
                on=CORNER_JOIN_KEYS
            )
        
        # Filter out boxes matched to themselves
        K_COUNT_DUMMY = "count"
        aligned = aligned[aligned[K_BOX_ID_FROM] != aligned[K_BOX_ID_TO]]
        aligned[K_COUNT_DUMMY] = 1

        # Each edge will appear twice but in reverse order. We ensure that
        # the edge with the lower identifier always appears in the 'from'
        # and the edge with the higher identifier in the 'to'
        # This operation will create duplicated rows.
        min_of_edges = aligned[[K_BOX_ID_FROM, K_BOX_ID_TO]].min(axis=1)
        max_of_edges = aligned[[K_BOX_ID_FROM, K_BOX_ID_TO]].max(axis=1)
        aligned[K_BOX_ID_FROM] = min_of_edges
        aligned[K_BOX_ID_TO] = max_of_edges

        # Two boxes (of the same depth) are connected only if they have 4 corners in common.
        # However, we created duplicated rows in the previous step, so the actual magic value
        # is 8 = 2*4
        edges = aligned.groupby(by=[K_BOX_ID_FROM, K_BOX_ID_TO, K_BOX_DEPTH])[K_COUNT_DUMMY].count()
        edges = edges[edges == 8]
        edges_df = edges.index.to_frame()
        min_depth = cast(int, edges_df[K_BOX_DEPTH].min())
        max_depth = cast(int, edges_df[K_BOX_DEPTH].max())

        depths = (corners.groupby(K_BOX_DEPTH)[K_BOX_SIZE].max() ** 3).to_dict()

        nodes_df = corners.groupby(K_BOX_ID)[[K_BOX_X, K_BOX_Y, K_BOX_Z, K_BOX_SIZE]].max()

        graphs = [
            FindCavitiesGraph(
                Graph.DataFrame(
                    edges_df[edges_df[K_BOX_DEPTH] == d].drop(K_BOX_DEPTH, axis=1),
                    use_vids=False,
                    directed=False
                ),
                d
            )
            for d in range(min_depth, max_depth + 1)
        ]

        return FindCavitiesResult(
            self,
            graphs,
            nodes_df,
            edges_df,
            depths
        )

    @staticmethod
    def find_cavities(
        points: NDArray[float64],
        empty_region_treshold: int
    ) -> FindCavitiesResult:
        ctx = FindCavitiesContext.construct(points, empty_region_treshold)
        corners = ctx.get_empty_corners()
        return ctx.construct_graph(corners)
    
def find_cavities(points: NDArray[float64], empty_region_treshold: int) -> FindCavitiesResult:
    return FindCavitiesContext.find_cavities(points, empty_region_treshold)

rot_360 = np.pi * 2
sphere_segs = 16

UNIT_SPHERE: NDArray[float64] = np.unique(
    np.concatenate([
        # Center of the sphere
        np.array([[0,0,0]]),

        # Points along the surface
        np.array([
            [
                np.sin(phi) * np.cos(theta),
                np.sin(phi) * np.sin(theta),
                np.cos(phi)
            ]
            for xa in range(sphere_segs - 1)
            for ya in range(sphere_segs - 1)
            for phi in [xa*rot_360/sphere_segs]
            for theta in [ya*rot_360/sphere_segs]
        ])
    ]),
    axis=0
)
    
def to_spheres_cloud(
    points: List[List[float]],
    radii: List[float],
    radii_scale: float
) -> NDArray[float64]:
    """Convert a list of points representing spheres ([cx,cy,cz,radius]) to a set
    of points positioned around the centers and surfaces of said spheres."""

    centers = np.array(points)
    radii_scales = np.array(radii) * radii_scale

    # Replicate the points of the sphere once per given radius scaling factor, then multiply
    # each of the points to its corresponding radius to scale the spheres
    spheres = np.repeat(UNIT_SPHERE, len(radii_scales), axis=0) * np.tile(radii_scales, (3, len(UNIT_SPHERE))).T
    #spheres = np.tile(UNIT_SPHERE, (len(radii_scales), 1)) * np.tile(radii_scales, (3, len(UNIT_SPHERE))).T

    # Replicate each center once per point in the sphere, then add them to the
    # points of the spheres in order to traslate them to the correct place
    return np.tile(centers, (len(UNIT_SPHERE),1)) + spheres