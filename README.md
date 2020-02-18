_________________________________________________________________________________

  ███╗   ██╗ ██████╗  ██████╗     ██████╗ ██╗      ██████╗ ██████╗  █████╗ ██╗         
  ████╗  ██║██╔═══██╗██╔════╝    ██╔════╝ ██║     ██╔═══██╗██╔══██╗██╔══██╗██║         
  ██╔██╗ ██║██║   ██║██║         ██║  ███╗██║     ██║   ██║██████╔╝███████║██║         
  ██║╚██╗██║██║   ██║██║         ██║   ██║██║     ██║   ██║██╔══██╗██╔══██║██║         
  ██║ ╚████║╚██████╔╝╚██████╗    ╚██████╔╝███████╗╚██████╔╝██████╔╝██║  ██║█████╗
  ╚═╝  ╚═══╝ ╚═════╝  ╚═════╝     ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚════╝    
                                                                                     
  ████████╗██╗██████╗ ███████╗    ███╗   ███╗ ██████╗ ██████╗ ███████╗██╗              
  ╚══██╔══╝██║██╔══██╗██╔════╝    ████╗ ████║██╔═══██╗██╔══██╗██╔════╝██║              
     ██║   ██║██║  ██║█████╗      ██╔████╔██║██║   ██║██║  ██║█████╗  ██║              
     ██║   ██║██║  ██║██╔══╝      ██║╚██╔╝██║██║   ██║██║  ██║██╔══╝  ██║              
     ██║   ██║██████╔╝███████╗    ██║ ╚═╝ ██║╚██████╔╝██████╔╝███████╗███████╗         
     ╚═╝   ╚═╝╚═════╝ ╚══════╝    ╚═╝     ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝╚══════╝         

01001110 01001111 01000011  01000111 01101100 01101111 01100010 01100001 01101100  
01010100 01101001 01100100 01100101  01001101 01101111 01100100 01100101 01101100 
_________________________________________________________________________________
|                                                                               |
|                        -~Version 1.0 (31/03/2020)~-                           |
|                      ~David Byrne (dbyrne@noc.ac.uk)~                         |
|_______________________________________________________________________________|

         o               .        ___---___                    .                   
                .              .--\        --.     .     .         .
                             ./.;_.\     __/~ \.     
                            /;  / `-'  __\    . \                            
         .         .       / ,--'     / .   .;   \        |
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

This repository contains the NEMO configuration information used for v1.0 of the NOC
global tide dataset. It contains the directories:

/ARCH/       - Architecture files used to run NEMO on ARCHER.

/EXP00/      - Experiment directory. Same as in EXP_R12

/EXP_R12/    - 1/12 degree experiment directory.

/EXP_R025/   - 1/4 degree experiment directory.

/MY_SRC/     - Modified NEMO source code files used for the configuration.

/WORK/       - NEMO source code used for the configuration.

/cpp_GTM.sh/ - compile options for NEMO configuration. 

The configuration is a modified version of the GO8p0.1 global NEMO configuration.
The Github wiki for this repository gives an overview of these changes. 
INPUT files are available and include:

domain_cfg.nc           -  domain file.

bt_tidal_dissipation.nc - tidal dissipation informations (open ocean).

ERA5_hourlymeans.nc     - Mean atmospheric forcing files.

Ask David Byrne (dbyrne@noc.ac.uk) for location of these files.
