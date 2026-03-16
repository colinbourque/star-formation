import numpy as np

with open(f'aux/wavelength_micron.inp','r') as f:
    wavs = f.readlines()
wavs = [float(wav) for wav in wavs[1:]]

aps = [2, 6, 20, 50]

with open('mod1step90/aperture_info.inp', 'w+') as f:
    f.write(f'1\n')
    f.write(f'100\n')
    for i in range(len(wavs)):
        if wavs[i] <= 10:
            f.write(f'{wavs[i]:.6e}\t{aps[0]:.6e}\n')
        elif 10 < wavs[i] <= 40:
            f.write(f'{wavs[i]:.6e}\t{aps[1]:.6e}\n')
        elif 40 < wavs[i] <= 100:
            f.write(f'{wavs[i]:.6e}\t{aps[2]:.6e}\n')
        elif 100 < wavs[i]:
            f.write(f'{wavs[i]:.6e}\t{aps[3]:.6e}\n')
