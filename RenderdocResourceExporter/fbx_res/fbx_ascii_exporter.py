'''
FBX ASCII 格式导出器
'''

import os
import sys
from pathlib import Path
from collections import defaultdict, OrderedDict
from textwrap import dedent

from PySide2 import QtWidgets, QtCore, QtGui

current_folder = Path(__file__).absolute().parent
father_folder = str(current_folder.parent)
sys.path.append(father_folder)

from qt_utilities.progress_bar_utility import ProgressBarUtility
from fbx_res.fbx_ascii_process_handler import FbxAsciiProcessHandler


class FbxAsciiExporter:

    @staticmethod
    def analyse_data(table_):

        model = table_.model()
        row_count = model.rowCount()
        column_count = model.columnCount()
        rows = range(row_count)
        columns = range(column_count)

        data = defaultdict(list)
        attr_list = set()

        for _, item in ProgressBarUtility.loop(columns, status_="收集mesh数据"):
            head = model.headerData(item, QtCore.Qt.Horizontal)
            values = [model.data(model.index(r, item)) for r in rows]
            if "." not in head:
                data[head] = values
            else:
                attr = head.split(".")[0]
                attr_list.add(attr)
                data[attr].append(values)

        for _, item in ProgressBarUtility.loop(attr_list, status_="重新排列mesh数据"):
            values_list = data[item]
            data[item] = [[float(values[r]) for values in values_list] for r in rows]

        return data, attr_list
    
    
    @staticmethod
    def export_fbx(save_path_, mapper_, table_):
        
        data, attr_list = FbxAsciiExporter.analyse_data(table_)

        POSITION = mapper_.get("POSITION")
        NORMAL = mapper_.get("NORMAL")
        BINORMAL = mapper_.get("BINORMAL")
        TANGENT = mapper_.get("TANGENT")
        COLOR = mapper_.get("COLOR")
        UV = mapper_.get("UV")
        UV2 = mapper_.get("UV2")
        ENGINE = mapper_.get("ENGINE")

        save_name = os.path.basename(os.path.splitext(save_path_)[0])
        ARGS = {
            "model_name": save_name,
            "LayerElementNormal": "",
            "LayerElementNormalInsert": "",
            "LayerElementBiNormal": "",
            "LayerElementBiNormalInsert": "",
            "LayerElementTangent": "",
            "LayerElementTangentInsert": "",
            "LayerElementColor": "",
            "LayerElementColorInsert": "",
            "LayerElementUV": "",
            "LayerElementUVInsert": "",
            "LayerElementUV2": "",
            "LayerElementUV2Insert": "",
        }

        # We'll decode the first three indices making up a triangle
        idx_dict = data["IDX"]
        value_dict = defaultdict(list)
        vertex_data = defaultdict(OrderedDict)

        for i, idx in enumerate(idx_dict):
            for attr in attr_list:
                value = data[attr][i]
                value_dict[attr].append(value)
                if idx not in vertex_data[attr]:
                    vertex_data[attr][idx] = value

        polygons = idx_dict
        if not polygons: return

        min_poly = min(polygons)
        idx_list = [str(idx - min_poly) for idx in idx_dict]
        idx_data = ",".join(idx_list)
        idx_len = len(idx_list)
        
        handler = FbxAsciiProcessHandler(
            {
                "POSITION": POSITION,
                "NORMAL": NORMAL,
                "BINORMAL": BINORMAL,
                "TANGENT": TANGENT,
                "COLOR": COLOR,
                "UV": UV,
                "UV2": UV2,
                "ENGINE": ENGINE,

                "ARGS": ARGS,

                "polygons": polygons,
                "min_poly": min_poly,
                "idx_list": idx_list,
                "idx_data": idx_data,
                "idx_len": idx_len,
                "idx_dict": idx_dict,
                "value_dict": value_dict,
                "vertex_data": vertex_data,
            })
        fbx = handler.run()

        with open(save_path_, "w") as f:
            f.write(dedent(fbx).strip())

