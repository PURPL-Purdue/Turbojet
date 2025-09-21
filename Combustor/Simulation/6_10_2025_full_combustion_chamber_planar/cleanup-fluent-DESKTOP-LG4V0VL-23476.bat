echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v252\fluent/ntbin/win64/winkill.exe"

start "tell.exe" /B "C:\PROGRA~1\ANSYSI~1\v252\fluent\ntbin\win64\tell.exe" DESKTOP-LG4V0VL 56161 CLEANUP_EXITING
timeout /t 1
"C:\PROGRA~1\ANSYSI~1\v252\fluent\ntbin\win64\kill.exe" tell.exe
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 19768) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 6924) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 22924) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 24276) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 23476) 
if /i "%LOCALHOST%"=="DESKTOP-LG4V0VL" (%KILL_CMD% 11984)
del "C:\Users\rjmcm\Documents\GitHub\Turbojet\Combustor\Simulation\6_10_2025_full_combustion_chamber_planar\cleanup-fluent-DESKTOP-LG4V0VL-23476.bat"
