'''
mesh 数据导出为 CSV 文件
'''

import os, csv, sys
from pathlib import Path

from PySide2 import QtWidgets, QtCore, QtGui

current_folder = Path(__file__).absolute().parent
father_folder = str(current_folder.parent)
sys.path.append(father_folder)

from qt_utilities.progress_bar_utility import ProgressBarUtility


class MeshToCsv:

    emgr = None

    def __init__(self, emgr_):
        self.emgr = emgr_
        pass

    def execute(self, table_, save_path_):

        # 先把旧的删掉
        if os.path.exists(save_path_):
            os.remove(save_path_)

        model = table_.model()
        row_count = model.rowCount()
        column_count = model.columnCount()
        rows = range(row_count)
        columns = range(column_count)

        head_names = list()
        row_data_list = list()

        for _, c_i in ProgressBarUtility.loop(columns, status_='Collect mesh data'):
            head = model.headerData(c_i, QtCore.Qt.Horizontal)
            row_data = [model.data(model.index(r_i, c_i)) for r_i in rows]

            head_names.append(head)
            row_data_list.append(row_data)

        with open(save_path_, "w", newline='') as csvFile:
            writer = csv.writer(csvFile)
            writer.writerow(head_names)

            for _, r_i in ProgressBarUtility.loop(rows, status_='Write to the csv file'):
                write_data = list()
                for row_data in row_data_list:
                    row_item = row_data[r_i]
                    write_data.append(row_item)
                writer.writerow(write_data)

        pass
