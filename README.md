    ███╗   ██╗ ██████╗  ██████╗     ██████╗ ██╗      ██████╗ ██████╗  █████╗ ██╗     
    ████╗  ██║██╔═══██╗██╔════╝    ██╔════╝ ██║     ██╔═══██╗██╔══██╗██╔══██╗██║     
    ██╔██╗ ██║██║   ██║██║         ██║  ███╗██║     ██║   ██║██████╔╝███████║██║     
    ██║╚██╗██║██║   ██║██║         ██║   ██║██║     ██║   ██║██╔══██╗██╔══██║██║     
    ██║ ╚████║╚██████╔╝╚██████╗    ╚██████╔╝███████╗╚██████╔╝██████╔╝██║  ██║███████╗
    ╚═╝  ╚═══╝ ╚═════╝  ╚═════╝     ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝                                                                            
      ████████╗██╗██████╗ ███████╗    ███╗   ███╗ ██████╗ ██████╗ ███████╗██╗          
      ╚══██╔══╝██║██╔══██╗██╔════╝    ████╗ ████║██╔═══██╗██╔══██╗██╔════╝██║          
         ██║   ██║██║  ██║█████╗      ██╔████╔██║██║   ██║██║  ██║█████╗  ██║          
         ██║   ██║██║  ██║██╔══╝      ██║╚██╔╝██║██║   ██║██║  ██║██╔══╝  ██║          
         ██║   ██║██████╔╝███████╗    ██║ ╚═╝ ██║╚██████╔╝██████╔╝███████╗███████╗     
         ╚═╝   ╚═╝╚═════╝ ╚══════╝    ╚═╝     ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝╚══════╝  
    _________________________________________________________________________________
    |                                                                               |
    |                        ~Version 1.0 (31/03/2020)~                           |
    |                      ~David Byrne (dbyrne@noc.ac.uk)~                         |
    |_______________________________________________________________________________|   
                .                                           o
                    .                                                      __
            o               .        ___---___                    .      / 0 (\          
                   .              .--\        --.     .     .         . |)  ~ o|
                                ./.;_.\     __/~ \.                      \ __ /
                               /;  /  -'  __\    . \                             
            .         .       / ,--'     / .   .;   \        |                    .
                             | .|       /       __   |      -O-       .
                            |__/    __ |  . ;   \ | . |      |
                            |      /  \\_    . ;| \___|    
               .    o       |      \  .~\\___,--'     |           .
                             |     | . ; ~~~~\_    __|
                |             \    \   .  .  ; \  /_/   .
               -O-        .    \   /         . |  ~/                  .
                |    .          ~\ \   .      /  /~          o
              .                   ~--___ ; ___--~       
                             .          ---         .              

This is the NEMO configuration used to generate Version 1.0 of the NOC global tide dataset.
The output from this NEMO configuration was used alongside data assimilation techniques 
to create a global dataset of tidal harmonic amplitudes and phases at 1/12 degree 
resolution. This repository contains the required NEMO files to install and run this 
configuration. The repository contains:

1. **/ARCH/**       - Architecture files used to run NEMO on ARCHER (now out of date).

2. **/EXP00/**      - Default experiment directory. This is a copy of the files in EXP_R12. 

3. **/EXP_R12/**    - Experiment directory for the 1/12-degree configuration.

4. **/EXP_R025/**   - Experiment directory for the 1/4-degree configuration.

5. **/WORK/**       - Base NEMO source code used for the configuration. This is based 
on the GO8p0.1 configuration and uses NEMO version 4.0.

6. **/MY_SRC/**     - Modified NEMO source code files used for this specific configuration. 
More information below.

7. **/cpp_GTM.sh/** - compile options for NEMO configuration. 

## Modifications in the MY_SRC directory:

**Overview.** Changes are to to do with updating the tidal forcing code slightly (according to [https://github.com/NOC-MSM/NEMO_cfgs/blob/master/recipes/docs/source/FES2014_NEMO.rst]), switching off surface fluxes and tracer processes as well the introduction as a basic internal wave drag parameterization from [reference]. Nico's changes are marked in the code with NB. My changes are marked with DB:.

By file:

1. tide.h90 :: Updated tidal potential forcing values (e.g. equilibrium tide amplitude). Data for more harmonic constituents than base NEMO. Should match the FES2014 dataset.
2. tide_mod.F90 :: Updated nodal factor equation for some harmonics.
3. sbctide.F90 :: Tidal potential equation for long period tides (previously not included).
4. dynspg_ts.F90 :: Reads tidal dissipation input array and applies to bottom friction component.
5. step.F90 :: Applies constant density by switching off calls to e.o.s, buoyancy frequency, tracer processes.
