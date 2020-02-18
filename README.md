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
