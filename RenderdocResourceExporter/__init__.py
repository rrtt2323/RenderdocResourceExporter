'''
资源导出插件的入口类
'''

import os

import qrenderdoc as qrd
from PySide2 import QtWidgets, QtCore, QtGui

from .qt_utilities.root_const_mapper import RootConstMapper as RCM

from .fbx_res.fbx_export_option_dialog import FbxExportOptionDialog
from .fbx_res.mesh_to_csv import MeshToCsv
from .fbx_res.csv_to_fbx import CsvToFbx


# 注册入口
def register(version_, pyrenderdoc_):
    emgr = pyrenderdoc_.Extensions()

    emgr.MessageDialog(RCM.get_text(RCM.c_welcome_1), RCM.get_text(RCM.c_welcome_2))

    if pyrenderdoc_.HasMeshPreview():
        emgr.RegisterPanelMenu(qrd.PanelMenu.MeshPreview, [RCM.get_text(RCM.c_export_fbx)], ExportFbx)


# 异常捕获
def error_log(func_):
    def wrapper(pyrenderdoc, data):
        emgr = pyrenderdoc.Extensions()
        try:
            func_(pyrenderdoc, data)
        except:
            import traceback
            exc = traceback.format_exc()
            emgr.ErrorDialog(exc, "Error!!!")
    
    return wrapper

# 导出fbx
@error_log
def ExportFbx(pyrenderdoc_, data_):
    emgr = pyrenderdoc_.Extensions()

    dialog = FbxExportOptionDialog(emgr) # 进行导出配置
    if not dialog.mqt.ShowWidgetAsDialog(dialog.init_ui()): return # 是不是主动取消的

    # 选择保存位置
    fbx_save_path = emgr.SaveFileName(RCM.get_text(RCM.c_save_fbx_file_to), '', '*.fbx')
    if not fbx_save_path:
        emgr.ErrorDialog(RCM.get_text(RCM.c_save_path_error), 'Error!!!')
        return
    # csv的地址
    csv_save_path = fbx_save_path.replace('.fbx', '.csv')

    # 直接从 QTableView 界面获取数据
    main_window = pyrenderdoc_.GetMainWindow().Widget()
    table = main_window.findChild(QtWidgets.QTableView, 'vsinData')
    model = table.model()
    row_count = model.rowCount()
    column_count = model.columnCount()
    rows = range(row_count)
    columns = range(column_count)
    if len(rows) <=1 and len(columns) <= 2:
        emgr.ErrorDialog(RCM.get_text(RCM.c_not_mesh_data), 'Error!!!')
        return

    # 导出CSV数据表
    meshToCsv = MeshToCsv(emgr)
    pyrenderdoc_.Replay().BlockInvoke(meshToCsv.execute(table, csv_save_path))
    # 读取CSV文件，导出FBX文件
    if os.path.exists(csv_save_path):
        csvToFbx = CsvToFbx(emgr)
        pyrenderdoc_.Replay().BlockInvoke(csvToFbx.execute(csv_save_path, fbx_save_path, dialog.export_config))

    emgr.MessageDialog(RCM.get_text(RCM.c_export_complete), '\(^o^)/')

    # 打开保存位置
    fbx_dir = os.path.dirname(fbx_save_path)
    if os.path.exists(fbx_dir):
        os.startfile(fbx_dir)

