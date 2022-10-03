'''
CSV 文件导出成 FBX 文件
'''

import os
import sys
from pathlib import Path

current_folder = Path(__file__).absolute().parent
father_folder = str(current_folder.parent)
sys.path.append(father_folder)

from qt_utilities.setting_utility import SettingUtility
from fbx_res.fbx_export_option_dialog_const_mapper import FbxExportOptionDialogConstMapper as FCM


class CsvToFbx:

    emgr = None

    def __init__(self, emgr_):
        self.emgr = emgr_
        pass

    def get_arg(self, config_, key_):
        if config_.get(key_):
            return '1'
        else:
            return '0'
        pass

    def execute(self, csv_path_, fbx_path_, export_config_):
        
        # 先把旧的删掉
        if os.path.exists(fbx_path_):
            os.remove(fbx_path_)

        # 读取设置
        export_normal = self.get_arg(export_config_, FCM.c_export_normal)
        export_tangent = self.get_arg(export_config_, FCM.c_export_tangent)
        export_uv = self.get_arg(export_config_, FCM.c_export_uv)
        export_uv2 = self.get_arg(export_config_, FCM.c_export_uv2)
        export_uv3 = self.get_arg(export_config_, FCM.c_export_uv3)
        
        # 拼接bat路径
        current_folder = Path(__file__).absolute()
        bat_path = str(current_folder.parent) + '\\csv_to_fbx.bat'

        # 运行参数格式化
        args = "%s %s %s %s %s %s %s " % (bat_path, csv_path_, export_normal, export_tangent, export_uv, export_uv2, export_uv3)
        # 通过执行命令行脚本来执行exe文件
        result = os.system(args)

        if result == 1:
            self.emgr.ErrorDialog('csv_to_fbx.bat return 1')

        pass

