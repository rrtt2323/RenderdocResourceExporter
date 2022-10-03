'''
FBX 导出选项设置窗 - 常量管理器
'''

import sys
from pathlib import Path

current_folder = Path(__file__).absolute().parent
father_folder = str(current_folder.parent)
sys.path.append(father_folder)

from qt_utilities.root_const_mapper import RootConstMapper
from qt_utilities.setting_utility import SettingUtility


class FbxExportOptionDialogConstMapper:

    # key值
    c_dialog_title = 'dialog_title'
    c_language = 'language'
    c_cancel = 'cancel'
    c_ok = 'ok'

    c_export_normal = 'export_normal'
    c_export_tangent = 'export_tangent'
    c_export_uv = 'export_uv'
    c_export_uv2 = 'export_uv2'
    c_export_uv3 = 'export_uv3'

    # 文本字典
    text_en = {
        c_dialog_title : 'FBX Export Option',
        c_language : 'Language',
        c_cancel : 'Cancel',
        c_ok : 'Ok',

        c_export_normal : 'Export Normal',
        c_export_tangent : 'Export Tangent',
        c_export_uv : 'Export UV',
        c_export_uv2 : 'Export UV2',
        c_export_uv3 : 'Export UV3',
    }
    text_cn = {
        c_dialog_title : 'FBX导出选项',
        c_language : '语言',
        c_cancel : '取消',
        c_ok : '确定',

        c_export_normal : '导出法线',
        c_export_tangent : '导出切线',
        c_export_uv : '导出 UV',
        c_export_uv2 : '导出 UV2',
        c_export_uv3 : '导出 UV3',
    }

    @staticmethod
    def get_text(key_):
        path, settings = SettingUtility.get_tempdir_setting(RootConstMapper.setting_preference)

        language = settings.value(RootConstMapper.c_language)

        if language == RootConstMapper.c_cn:
            return FbxExportOptionDialogConstMapper.text_cn[key_]
        else:
            return FbxExportOptionDialogConstMapper.text_en[key_]
        
        pass
