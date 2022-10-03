'''
进度条实用工具
'''

from PySide2 import QtWidgets, QtCore, QtGui


class ProgressBarUtility(QtWidgets.QProgressDialog):
    
    #构造函数
    def __init__(self,
        status_ = u"progress...",
        button_text_ = u"Cancel",
        minimum_ = 0,
        maximum_ = 100,
        parent_ = None,
        title_ = "",
    ):
        super(ProgressBarUtility, self).__init__(parent_)

        self.setWindowFlags(self.windowFlags())
        self.setWindowModality(QtCore.Qt.WindowModal)
        self.setWindowTitle(status_ if title_ else title_)

        bar = QtWidgets.QProgressBar(self)
        bar.setStyleSheet(
            """
            QProgressBar {
                color:white;
                border: 1px solid black;
                background: gray;
            }

            QProgressBar::chunk {
                background: QLinearGradient( x1: 0, y1: 0, x2: 1, y2: 0,
                stop: 0 #78d,
                stop: 0.4999 #46a,
                stop: 0.5 #45a,
                stop: 1 #238 );
                border: 1px solid black;
            }
            """
        )
        bar.setAlignment(QtCore.Qt.AlignCenter)

        self.setBar(bar)
        self.setLabelText(status_)
        self.setCancelButtonText(button_text_)
        self.setRange(minimum_, maximum_)
        self.setValue(minimum_)
        
        # NOTE show the progressbar without blocking
        self.show()
        QtWidgets.QApplication.processEvents()


    #补充构造函数
    @classmethod
    def loop(cls, seq_, **kwargs_):
        self = cls(**kwargs_)
        
        #如果参数里没有规定最大值，就用序列的长度
        if not kwargs_.get("maximum"):
            self.setMaximum(len(seq_))

        for i, item in enumerate(seq_, 1):
            if self.wasCanceled():
                break
            try:
                yield i, item # with body executes here
            except:
                import traceback
                traceback.print_exc()
                self.deleteLater()
            self.setValue(i)
        
        self.deleteLater()
