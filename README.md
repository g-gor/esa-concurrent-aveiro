# esa-concurrent-aveiro

main.m executes a simple Lambert Battin algorithm to find the best date to change orbit from Earth to the comet 311P, minimizing delta-v.
Another comets can be introduced by adding its position and velocity on the initial date 2032/12/31

main_porkchops.m executes the Lambert Battin algorithm to display visually the best timeframes to do a orbit transfer, because this plot records the required delta-v as a function of time.
This script sweeps trough a 6 year launch timeframe and also a 6 year arrival timeframe. The arrival date is at least 1 day after the launch date, which does not make sense for long trips, but this can be adjusted in the script by changing the lr_offset
