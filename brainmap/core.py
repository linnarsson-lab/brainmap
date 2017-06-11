import numpy as np
import logging
from typing import *
from allensdk.api.queries.ontologies_api import OntologiesApi
import matplotlib.pyplot as plt
from colormap import Color
try:
    from ipywidgets import interact, interactive, fixed, interact_manual
    import ipywidgets as widgets
except ImportError:
    logging.warn("ipywidgets is not installed some function will not be available")


class AllenBrainStructure:
    def __init__(self, object_dict: Dict[str, Any], atlas: Any) -> None:
        for k, v in object_dict.items():
            setattr(self, k, v)
        self.children = []
        self.parent = None
        self._atlas = atlas
        
    def add_children(self, other: Any) -> None:
        if other not in self.children:
            self.children.append(other)
            other.set_parent(self)
            
    def set_parent(self, other):
        if self.parent is None:
            self.parent = other
            other.add_children(self)
            
    def connect(self):
        if self.parent_structure_id is not None:
            parent_obj = self._atlas[self.parent_structure_id]
            self.set_parent(parent_obj)
        else:
            pass
    
    def __repr__(self):
        head = super(AllenBrainStructure, self).__repr__()
        memory_address = head.split(" ")[-1][:-1]
        return '<AllenBrainStructure %s "%s" at %s>' % (self.id, self.safe_name, memory_address)
    
    def _repr_html_(self):
        html = "<p>"
        html += "<table>"
        html += '<tr align="center">'
        html += "<td>&nbsp;</td>"
        html += "<td><strong>%s</strong></td>" % self.safe_name
        html += "</tr>"
        entries_to_print = ["acronym", "name", "id", "graph_order",
                            "st_level", "depth", "structure_id_path"]
        for k in entries_to_print:
            v = getattr(self, k)
            html += '<tr align="center">'
            html += "<td><strong>" + str(k) + "</strong></td>"
            html += "<td>" + str(v) + "</td>"
            html += "</tr>"
        html += '<tr><td><strong>color</strong></td>'
        html += '<td bgcolor="#%s">%s</td></tr>' % (self.color_hex_triplet, self.color_hex_triplet)
        html += "</table>"
        return html


class AllenBrainReference:
    def __init__(self, graph="adult"):
        self._onto = OntologiesApi()
        if graph == "adult":
            self._struct_dicts_list = self._onto.get_structures(structure_graph_ids=1, num_rows='all')
        elif graph == "development":
            self._struct_dicts_list = self._onto.get_structures(structure_graph_ids=17, num_rows='all')
        self._structures = {i["id"]: AllenBrainStructure(i, self) for i in self._struct_dicts_list}
        self._link()
    
    def __getitem__(self, key: int) -> Any:
        return self._structures[key]
    
    def __iter__(self) -> Any:
        for k, v in self._structures.items():
            yield v
    
    def _link(self) -> None:
        for name, structure in self._structures.items():
            structure.connect()


class Reference3D:
    def __init__(self):
        pass


class ColoredVolumetric:
    def __init__(self, allen_vol_data):
        self.vol_data = allen_vol_data
        # Color(i) for self.vol_data.color_list
        
    def __getitem__(some_slice):
        sliced = self.vol_data[some_slice]
        sliced = np.tile(sliced[:, :, None], 3)
        
    
class AllenVolumetricData:
    def __init__(self, filename: str, reference: Any=None) -> None:
        self.filename = filename
        self.zip_container = zipfile.ZipFile(filename)
        for i in zip_container.infolist():
            if ".mhd" in i.filename:
                mhd_file = i
            if '.raw' in i.filename:
                raw_file = i
        info_file = self.zip_container.open(mhd_file).read().decode("ascii")
        entries_file = [i.split(" = ") for i in info_file.rstrip().split("\n")]
        self.file_info = {k: ([int(i) for i in v.split(" ")] if (" " in v) else (k, v)) for k, v in entries_file}
        self.shape = tuple(file_info['DimSize'])
        self.is_label = ("UINT" in file_info['ElementType'][1])
        self.file_type = 'uint32' if self.is_label else 'uint8'
        logging.debug("Reading data file")
        buffer = zip_container.open(raw_file).read()
        array1d = np.fromstring(buffer, dtype=file_type)
        if self.is_label:
            self.ids, array1d = np.unique(array1d, return_inverse=True)
        self._values = array1d.reshape(dim_size, order='F')
        self.zip_container.close()
        self.reference = reference
        if self.is_label:
            self._colored = ColoredVolumetric(self)
        
    def __getitem__(self, slice_obj):
        return self._values[slice_obj]
    
    @property
    def reference(self):
        try:
            return self._reference
        except NameError:
            self._reference = None
            return self.reference
    
    @reference.setter
    def reference(self, reference_object):
        self._reference = reference_object
        self.color_list = []
        for idx in self.ids:
            try:
                self.color_list.append("#" + self._reference[idx].color_hex_triplet)
            except KeyError:
                self.color_list.append("#000000")
    
    @property    
    def colored(self):
        return self._colored
    
    def plot_slides(self, coronal, sagittal, ss=None, fig=None, return_figure=False):
        if fig is None:
            fig = plt.gcf()
        if ss is None:
            ss = plt.GridSpec(1, 1)[0]
        if return_figure:
            fig.clear()
        if self.is_label and self.reference:
            _cmap = ListedColormap(self.color_list)
        else:
            _cmap = "gray"
        gs = GridSpecFromSubplotSpec(1, 2, subplot_spec=ss)
        ax = fig.add_subplot(gs[0])
        ax.imshow(self[coronal, :, :], cmap=_cmap)
        ax.axvline(sagittal, c='y', lw=1)
        ax = fig.add_subplot(gs[1])
        ax.imshow(self[:, :, sagittal].T, cmap=_cmap)
        ax.axvline(coronal, c='y', lw=1)
        if return_figure:
            return fig
        
    def interactive_slides(self):
        fig = plt.figure(figsize=(12, 5))
        gs = plt.GridSpec(1, 1)
        interact(plot_slides,
                 coronal=(0, array_3d.shape[0] - 1),
                 sagittal=(0, array_3d.shape[-1] - 1),
                 ss=fixed(gs[0]), fig=fixed(fig))
        plt.clf()
