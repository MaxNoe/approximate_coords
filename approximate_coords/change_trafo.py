from contextlib import contextmanager
from astropy.coordinates import frame_transform_graph


@contextmanager
def use_transformation(fromsys, tosys, transform):
    old_transform = frame_transform_graph._graph[fromsys].get(tosys)

    try:
        frame_transform_graph.add_transform(fromsys, tosys, transform)
        yield
    finally:
        frame_transform_graph.remove_transform(fromsys, tosys, transform)
        if old_transform is not None:
            frame_transform_graph.add_transform(fromsys, tosys, old_transform)
