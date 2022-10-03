'''
根 - 常量管理器
'''

import sys
from pathlib import Path

# 切换到上一层目录，以便导入父目录的模块
current_folder = Path(__file__).absolute().parent
father_folder = str(current_folder.parent)
sys.path.append(father_folder)

from qt_utilities.setting_utility import SettingUtility


class RootConstMapper:

    setting_preference = 'renderdoc_resource_exporter_preference.ini'

    c_language = 'language'
    c_en = 'en'
    c_cn = 'cn'

    c_welcome_1 = 'welcome_1'
    c_welcome_2 = 'welcome_2'
    c_export_fbx = 'export_fbx'
    c_save_fbx_file_to = 'save_fbx_file_to'
    c_save_path_error = 'save_path_error'
    c_not_mesh_data = 'not_mesh_data'
    c_export_complete = 'export_complete'

    text_en = {
        c_welcome_1 : 'The resource exporter was registered successfully.',
        c_welcome_2 : 'Hey! Cowboy!',
        c_export_fbx : 'Export FBX',
        c_save_fbx_file_to : 'Save FBX File To',
        c_save_path_error : 'Save Path Error!!!',
        c_not_mesh_data : 'Not Mesh Data!!!',
        c_export_complete : 'Export complete. The save location then opens automatically.',
    }
    text_cn = {
        c_welcome_1 : '资源导出器注册成功',
        c_welcome_2 : '嘿！牛仔！',
        c_export_fbx : '导出FBX',
        c_save_fbx_file_to : '保存FBX文件到',
        c_save_path_error : '保存路径错误!!!',
        c_not_mesh_data : '没有网格数据!!!',
        c_export_complete : '导出完成。自动打开保存位置。',
    }

    @staticmethod
    def get_text(key_):
        path, settings = SettingUtility.get_tempdir_setting(RootConstMapper.setting_preference)

        language = settings.value(RootConstMapper.c_language)

        if language == RootConstMapper.c_cn:
            return RootConstMapper.text_cn[key_]
        else:
            return RootConstMapper.text_en[key_]
        
        pass
