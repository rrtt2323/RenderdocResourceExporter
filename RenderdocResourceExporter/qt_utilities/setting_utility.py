'''
配置文件实用工具
'''

import os, tempfile

from PySide2 import QtWidgets, QtCore, QtGui

class SettingUtility:

    @staticmethod
    def get_tempdir_setting(iniName):
        path = os.path.join(tempfile.gettempdir(), iniName)
        settings = QtCore.QSettings(path, QtCore.QSettings.IniFormat)
        return path, settings
        pass
