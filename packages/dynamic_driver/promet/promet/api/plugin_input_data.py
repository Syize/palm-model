import pyproj


class PluginInputData:

    attributes = dict()
    available_time_slots = list()

    def load_grids(
            self,
            transformer: pyproj.Transformer,
            bounds_x: tuple,
            bounds_y: tuple,
            bounds_z: tuple,
            bounds_zsoil: tuple,
    ):
        pass

    def discover_dataset(self):
        pass

    def get_variable_timeframe(
            self,
            input_variable,
            time_index,
            terrain_extended=False,
    ):
        pass

    def get_terrain_frame(self, time_index):
        pass

    def get_metadata(self) -> dict:
        pass
