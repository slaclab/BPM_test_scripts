clear; clf;
Data_directory = '/afs/slac/g/lcls/users/BPM/LCLS_II/Data/';
dir_name  = datestr(now,'mm_dd_yyyy_HHMMSS');
directory = Data_directory;

concate_filename = char(strcat('BPM_Res_data_SNR_168_',dir_name,'.mat'));
edef = '30';
tpg_pv = strcat('BSA:B084:2:', edef, ':');

% verify calibration triggers are disabled
lcaPut('BPMS:B084:200:CALBTCTL', 0);

lcaClear
mat_fname = [Data_directory,concate_filename];
nshots = input('How many sets of data was acquired? ');
%nshots=2000;
lcaPut(strcat(tpg_pv, 'MEASCNT'), nshots);

attn1_val=6;
lcaPut('BPMS:B084:200:ATT1', attn1_val);

attn2_val=10;
lcaPut('BPMS:B084:200:ATT2', attn2_val);

pause(1);
raw_wave=lcaGet('BPMS:B084:200:RWAV');
pause(.5);
lcaPut(strcat(tpg_pv, 'RATEMODE'), 0);
pause(.5)
lcaPut(strcat(tpg_pv, 'FIXEDRATE'), 1)
pause(.5);
lcaPut(strcat(tpg_pv, 'MEASSEVR.VAL'), 2);
pause(0.5);
lcaPut(strcat(tpg_pv, 'DESTMODE'), 0);
pause(0.5);

pause(.5);
lcaPut(strcat(tpg_pv, 'CTRL'), 1);
pause(.5);
p = lcaGet(strcat(tpg_pv, 'CTRL.RVAL'));
while p==1
  p = lcaGet(strcat(tpg_pv, 'CTRL.RVAL'));
  pause(1)
end
y_data=lcaGet(strcat('BPMS:B084:200:YHST', edef));
x_data=lcaGet(strcat('BPMS:B084:200:XHST', edef));
res_y=std(y_data(1:nshots));
res_x=std(x_data(1:nshots));
mat_fname
yres=res_y*1e3
xres=res_x*1e3
figure(1)
plot(x_data(1:nshots), y_data(1:nshots), '*')
title('Position of beam');
xlabel('X position (um)')
ylabel('Y position (um)')
figure(2)
max_red=max(raw_wave(1:127));
max_yellow=max(raw_wave(128:255));
max_blue=max(raw_wave(256:384));
max_green=max(raw_wave(384:512));
plot(raw_wave)

save(mat_fname);
