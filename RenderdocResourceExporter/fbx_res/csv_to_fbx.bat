:: 用来调用 RenderdocCSVToFBX.exe 的批处理文件

::cd /d %~dp0
set current_path=%~dp0

echo current_path=%current_path%
start %current_path%RenderdocCSVToFBX.exe %1 %2 %3 %4 %5 %6

::pause
