PROGRAM DigiBog_Hydro

!Version: 17/10/2024

!== Software licence ===========================================================
!DigiBog_Hydro. Copyright (c) 2013  See authors details below.

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the gnu General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but without any warranty; without even the implied warranty of
!    merchantability or fitness for a particular purpose.  See the
!    gnu General Public License for more details.

!    You should have received a copy of the gnu General Public License
!    along with this program.  If not, see https://www.gnu.org/licenses/
!===============================================================================

!-------------------------------------------------------------------------------
! Section 1.0 Program header
!-------------------------------------------------------------------------------

  !Description
  !A model to simulate water tables in ombrotrophic bogs and other shallow
  !aquifers using the Boussinesq equation. The model uses centimetres and
  !seconds. Details of how the model works may be obtained from the code
  !owners (email addresses below).

  !Current code owners
  !Paul J. Morris, Andy J. Baird*, and Lisa R. Belyea
  !*School of Geography
  !University of Leeds
  !Leeds
  !LS2 9JT
  !p.j.morris@reading.ac.uk
  !a.j.baird@leeds.ac.uk
  !l.belyea@qmul.ac.uk

  !Modification history of code
  !Programmer           Date           Modifications
  !============================================================================
  !Andrew J. Baird      08/04/2005     Original 1.5-d code
  !----------------------------------------------------------------------------
  !Paul J. Morris       24/02/2006     Conversion to 2.5-d and Fortran
  !                                    Standards implemented
  !----------------------------------------------------------------------------
  !Paul J. Morris       27/02/2006     Testing, minor corrections
  !----------------------------------------------------------------------------
  !Paul J. Morris       05/03/2006     Subroutines 'column_activation' and
  !                                    'steady_state_check' written. Former
  !                                    no longer exists (removed 18/07/2013).
  !----------------------------------------------------------------------------
  !Paul J. Morris       19/06/2006     Testing completed
  !----------------------------------------------------------------------------
  !Paul J. Morris       20/03/2007     Code cleaning
  !----------------------------------------------------------------------------
  !Paul J. Morris       24/04/2007     Further cleaning, including removal of
  !                                    replicated spatial step reference in
  !                                    'move_water' subroutine
  !----------------------------------------------------------------------------
  !Paul J. Morris       09/05/2007     Above-ground storage facilitated in
  !                                    2.5-d version
  !----------------------------------------------------------------------------
  !Paul J. Morris       02/07/2008     Final code cleaning, full annotation
  !----------------------------------------------------------------------------
  !Andy J. Baird        18/07/2013     Addition of variable net rainfall read
  !                                    in from file. Change to how boundary
  !                                    conditions specified; zero-flow
  !                                    Neumann condition also now allowed for,
  !                                    as are internal boundaries. Change to
  !                                    how output written to file (can now
  !                                    write to file multiple times before end
  !                                    of run).
  !----------------------------------------------------------------------------
  !Andy J. Baird        27/08/2013     Code cleaning and de-bugging of version
  !                                    from 18/07/2013
  !----------------------------------------------------------------------------
  !Dylan M. Young       17/10/2024     Added aet function to hydro_procedures
  !----------------------------------------------------------------------------

  !Model data

  !Input files
  !All input files are read by the main program, before any subroutines are
  !called. These input files must be in the same folder as the project file
  !and all source files. The input files are:

  !hydro_parameters.txt
  !This file contains all the scalar parameters used in the model

  !hydro_no_layers.txt
  !This file contains the number of layers to be used in each vertical,
  !square-sectioned column. Sorted in order of x, then y (i.e., for each x
  !poaition, values for all y positions are given before moving onto next x).

  !hydro_thickness.txt
  !This file contains the vertical thickness of each peat layer which makes
  !up the model aquifer. Sorted in order of x, then y, then z (i.e., one
  !layer at a time, from bottom to top of each column).

  !hydro_k.txt
  !This file contains the hydraulic conductivity of each peat layer which
  !makes up the model aquifer. Sorted in order of x, then y, then z (i.e., one
  !layer at a time, from bottom to top of each column).

  !hydro_s.txt
  !This file contains the drainable porosity of each peat layer which makes
  !up the model aquifer. Sorted in order of x, then y, then z (i.e., one
  !layer at a time, from bottom to top of each column).

  !hydro_column_status.txt
  !This files contains information on whether a column or cell is active
  !("on"), inactive ("off"), a Dirichlet boundary ("diri"), or a Neumann
  !boundary ("neu")

  !hydro_wt_bc.txt
  !This file contains the values for initial water-table height or the
  !boundary condition in each column. Sorted in order of x, then y.

  !hydro_baltitude.txt
  !This file contains the base altitude for each column, allowing the user
  !to represent a sloping or uneven substrate if necessary. Sorted in
  !order of x, then y.

  !hydro_rainfall.txt
  !This file contains daily rainfall for each day of the simulation

  !Output files
  !Three output files are used, as follows

  !hydro_wt_output.txt
  !This file records the water-table position, relative to the impermeable
  !base, at the end of each write-out period

  !hydro_x_indices.txt
  !This file records the x position of each column, at the end of the
  !model run, in order to allow the user to plot the water-table data in an
  !(x, y, z) triplet

  !hydro_y_indices.txt
  !This file records the y position of each column, at the end of the
  !model run, in order to allow the user to plot the water-table data in an
  !(x, y, z) triplet

  !Code description
  !Language: Fortran 95
  !Software standards: "UK Meteorological Office Fortran 90 Standards".
  !http://research.metoffice.gov.uk/research/nwp/numerical/fortran90/
          !f90_standards.html

  !Declarations

  !Modules used
  USE hydro_procedures

  IMPLICIT NONE

  !Global/local scalars
  INTEGER :: alloc_error, &     !Internal array allocation error flag
             t_extent, &        !Length of simulation in timesteps
             output_interval, & !Output interval in timesteps
             daily_timesteps, & !Number of timesteps per day
             time_counter, &    !Counter for use in main model loop
             output_counter, &  !Output counter for main model loop
             rain_counter, &    !Counter to check for update to rainfall
             x_extent, &        !Model grid extent in x direction
             y_extent, &        !Model grid extent in y direction
             z_extent, &        !Max number of layers per column
             x, &               !Spatial index counter
             y, &               !Spatial index counter
             z, &               !Spatial index counter
             steady_columns     !Threshold number of steady columns

  REAL(KIND=q) :: timestep, &         !Model temporal hincrement
                  spatial_step, &     !Model horizontal increment
                  rainfall, &         !Rainfall
                  pet,      &         !Potential evapotranspiration
                  n_aet,    &         !Defines the shape of aet extinction fun
                  aet_extinct, &      !Depth at which aet = 0
                  pond_depth, &       !Depth of surface ponding
                  steady_threshold    !Steady-state criterion per column

  CHARACTER(LEN=1) :: param_error  !Internal parameter error flag

  CHARACTER(LEN=9) :: hydro_status !Indicates model is steady or transient

  CHARACTER(LEN=25) :: data_file_name_1, &  !Parameter input file
                       data_file_name_2, &  !No. layers input file
                       data_file_name_3, &  !Layer thickness input file
                       data_file_name_4, &  !Layer K input file
                       data_file_name_5, &  !Layer s input file
                       data_file_name_6, &  !Activation status file
                       data_file_name_7, &  !Column WTs / BCs input file
                       data_file_name_8, &  !Base altitude input file
                       data_file_name_9, &  !Rainfall rate for each day
                       data_file_name_10, & !WT output file
                       data_file_name_11, & !x coordinates output file
                       data_file_name_12, & !y coordinates output file
                       data_file_name_13, & !Potential evapotranspiration
                       data_file_name_14    !AET status for gully water

  !Global/local arrays
  logical, allocatable, dimension(:, :) :: is_gully

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: no_layers !Number of layers per column

  REAL(KIND=q), ALLOCATABLE, DIMENSION(:,:) :: base_altitude, & !Above datum
                                               water_change, &
                                               water_table, &   !Above base
                                               wk_mean          !Depth-av. K

  !layer_attributes stores layer thickness, K and s
  !transmissivity stores layer elevation above base and transmissivity
  !layer_storage stores layer elevation above base and each layer's water
  !capacity as volume per unit area; i.e. expressed as a depth
  REAL(KIND=q), ALLOCATABLE, DIMENSION(:,:,:,:) :: layer_attributes, &
                                                   transmissivity, &
                                                   layer_storage

  !activation_status indicates whether/how column participates in simulation
  CHARACTER(8), ALLOCATABLE, DIMENSION(:,:) :: activation_status
  integer, ALLOCATABLE, DIMENSION(:,:) :: gully_status


!-------------------------------------------------------------------------------
! Section 2.0 Data input and initialisation
!-------------------------------------------------------------------------------

  !Name input data files and output data file
  data_file_name_1  = "hydro_parameters.txt"
  data_file_name_2  = "hydro_no_layers.txt"
  data_file_name_3  = "hydro_thickness.txt"
  data_file_name_4  = "hydro_k.txt"
  data_file_name_5  = "hydro_s.txt"
  data_file_name_6  = "hydro_column_status.txt"
  data_file_name_7  = "hydro_wt_bc_input.txt"
  data_file_name_8  = "hydro_baltitude.txt"
  data_file_name_9  = "hydro_rainfall.txt"
  data_file_name_10 = "hydro_wt_output.txt"
  data_file_name_11 = "hydro_x_indices.txt"
  data_file_name_12 = "hydro_y_indices.txt"
  data_file_name_13 = "hydro_pet.txt"
  data_file_name_14 = "hydro_aet_status.txt"

  !Open files for input of data and file for output of data
  open(unit=10,  file=data_file_name_1,  status="old")
  open(unit=20,  file=data_file_name_2,  status="old")
  open(unit=30,  file=data_file_name_3,  status="old")
  open(unit=40,  file=data_file_name_4,  status="old")
  open(unit=50,  file=data_file_name_5,  status="old")
  open(unit=60,  file=data_file_name_6,  status="old")
  open(unit=65,  file=data_file_name_14, status="old")
  open(unit=70,  file=data_file_name_7,  status="old")
  open(unit=80,  file=data_file_name_8,  status="old")
  open(unit=90,  file=data_file_name_9,  status="old")
  open(unit=100, file=data_file_name_10, status="replace")
  open(unit=110, file=data_file_name_11, status="replace")
  open(unit=120, file=data_file_name_12, status="replace")
  open(unit=130, file=data_file_name_13, status="old")


  !Licence statement -----------------------------------------------------------
  print *
  print *,  "DigiBog_Hydro  Copyright (c) 2013 the authors. ", &
    "This program comes with absolutely NO WARRANTY; for details see the ", &
    "included licence. This is free software, and you are welcome to ", &
    "redistribute it under certain conditions; for details see the ", &
    "included licence."
  print *
  !-----------------------------------------------------------------------------
  print *, "Press Enter to accept and continue"
  read( *, * )


  !Read data from parameter file
  READ (10, *) t_extent
  READ (10, *) output_interval
  READ (10, *) daily_timesteps
  READ (10, *) x_extent
  READ (10, *) y_extent
  READ (10, *) z_extent
  READ (10, *) steady_columns
  READ (10, *) spatial_step
  READ (10, *) pond_depth
  READ (10, *) steady_threshold
  read (10, *) n_aet
  read (10, *) aet_extinct

  !Error check on data
  WRITE (*, *) "The following parameter values have been read from ", &
                data_file_name_1
  WRITE (*, '(A19, I6, A7)')  "runtime =          ", t_extent, " tsteps"
  WRITE (*, '(A19, I4, A7)')  "output interval =  ", output_interval, " tsteps"
  WRITE (*, '(A19, I4, A7)')  "tsteps per day =   ", daily_timesteps, " tsteps"
  WRITE (*, '(A19, I4, A8)')  "x extent =         ", x_extent, " columns"
  WRITE (*, '(A19, I4, A8)')  "y extent =         ", y_extent, " columns"
  WRITE (*, '(A19, I4, A7)')  "z extent =         ", z_extent, " layers"
  WRITE (*, '(A19, I5)')      "steady columns =   ", steady_columns
  WRITE (*, '(A19, F7.2, A3)')"spatial step =     ", spatial_step, " cm"
  WRITE (*, '(A19, F7.2, A3)')"ponding depth =    ", pond_depth, " cm"
  WRITE (*, '(A19, F6.3, A3)')"steady threshold = ", steady_threshold, " cm"
  WRITE (*, '(A19, F4.2, A9)')"AET ext. param. =   ", n_aet, " unitless"
  WRITE (*, '(A19, f6.2, A3)')"AET ext. depth =    ", aet_extinct, " cm"
  WRITE (*, '(A31)') "Are these values correct (Y/N)?"
  READ *, param_error
  SELECT CASE (param_error)
  CASE ("n")
    WRITE (*, *) "Error in data - model terminating"
    STOP
  CASE ("N")
    WRITE (*, *) "Error in data - model terminating"
    STOP
  CASE DEFAULT
    WRITE (*, *) "Data appear to be okay - model run continues"
  END SELECT

  !Allocate length to arrays
  ALLOCATE(no_layers(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for no_layers"
    STOP
  END IF

  ALLOCATE(base_altitude(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for base_altitude"
    STOP
  END IF

  ALLOCATE(water_change(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for water_change"
    STOP
  END IF

  ALLOCATE(water_table(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for water_table"
    STOP
  END IF

  ALLOCATE(wk_mean(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for wk_mean"
    STOP
  END IF

  ALLOCATE(layer_attributes(x_extent, y_extent, z_extent, 3), &
                                                      STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for layer_attributes"
    STOP
  END IF

  ALLOCATE(transmissivity(x_extent, y_extent, z_extent, 2), STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for transmissivity"
    STOP
  END IF

  ALLOCATE(layer_storage(x_extent, y_extent, z_extent, 2), STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for layer_storage"
    STOP
  END IF

  ALLOCATE(activation_status(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for activation_status"
    STOP
  END IF

  ALLOCATE(gully_status(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for gully status"
    STOP
  END IF

  ALLOCATE(is_gully(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for is_gully"
    STOP
  END IF

  !Read data from files to rank-two arrays
  DO x = 1, x_extent
    DO y = 1, y_extent
      READ (20, *) no_layers(x, y)
      !Incorporate above-ground storage layer
      no_layers (x, y) = no_layers(x, y) + 1
      READ (60, *) activation_status(x, y)
      READ (65, *) gully_status(x, y)
      READ (70, *) water_table(x, y)
      READ (80, *) base_altitude(x, y)
    END DO
  END DO

  !Initialise layer_attributes
  layer_attributes = 0.0

  !Initialise is_gully
  is_gully = .false.
  is_gully = gully_status == 1

  DO x = 1, x_extent
    DO y = 1, y_extent
      !Initialise above-ground storage layer properties
      layer_attributes(x, y, no_layers(x, y), 1) = pond_depth
      layer_attributes(x, y, no_layers(x, y), 2) = 5
      layer_attributes(x, y, no_layers(x, y), 3) = 0.9
    END DO
  END DO

  !Read data from layer files to rank-four array
  DO x = 1, x_extent
    DO y = 1, y_extent
      DO z = 1, (no_layers(x, y) - 1)
        READ (30, *) layer_attributes(x, y, z, 1)
        READ (40, *) layer_attributes(x, y, z, 2)
        READ (50, *) layer_attributes(x, y, z, 3)
      END DO
    END DO
  END DO

  !Initialise remaining arrays
  water_change = 1.0
  wk_mean = 1.0
  transmissivity = 1.0
  layer_storage = 1.0

  !Initialise elapsed time
  time_counter = 0

  !Initialise check counter for write out
  output_counter = 0

  !Calculate value of timestep in seconds
  timestep = (1.0/(REAL(daily_timesteps))) * 86400.0

  !Read first net rainfall rate
  READ (090, *) rainfall
  READ (130, *) pet

  !Initialise rain counter
  rain_counter = 0


!-------------------------------------------------------------------------------
! Section 3.0 Calls to subroutines; time management; output
!-------------------------------------------------------------------------------
  !Calculate transmissivity profile for each column
  CALL trans_height(x_extent, y_extent, no_layers, layer_attributes, &
                    transmissivity, activation_status)

  !Calculate amount of water (expressed as a depth) which may be stored in each
  !layer within each column
  CALL layer_water_depth(x_extent, y_extent, no_layers, layer_attributes, &
                         layer_storage, activation_status)


  !Start of hydrological loop
  DO

    !(Re-)calculate the depth-averaged K below the water table for each column
    CALL wat_k_mean(x_extent, y_extent, no_layers, water_table, wk_mean, &
                    layer_attributes, transmissivity, activation_status)

    !Calculate the amount of water (expressed as a depth) that moves between
    !each column. The flow law is based on the Boussinesq equation.
    CALL move_water(x_extent, y_extent, timestep, spatial_step, rainfall, &
                        base_altitude, water_change, water_table, wk_mean, &
                        activation_status, no_layers, layer_storage, pet, &
                        n_aet, aet_extinct, is_gully)

    !Check for steady-state hydrological behaviour
    !CALL steady_state_check(x_extent, y_extent, steady_columns, &
                        !steady_threshold, hydro_status, water_change, &
                            !activation_status)

    !Update position of the water table in each column
    CALL water_table_update(x_extent, y_extent, no_layers, water_change, &
                            water_table, layer_attributes, layer_storage, &
                            activation_status)

    !Check for model write out
    output_counter = output_counter + 1
    IF (output_counter == output_interval) THEN
      !Re-set output counter
      output_counter = 0
      !Write results to file
      DO x = 1, x_extent
        DO y = 1, y_extent
          IF(activation_status (x,y) == "off" &
             .OR. activation_status (x,y) == "diri" &
             .OR. activation_status (x,y) == "neu") THEN
            WRITE(100, '(20F20.8)') -999.0
          ELSE
            WRITE(100, '(20F20.8)') water_table(x,y)
          END IF
        END DO
      END DO
    END IF

    !Update elapsed time and check for model termination
    time_counter = time_counter + 1
    IF (time_counter == t_extent) THEN
      EXIT
    END IF

    !Update rainfall rate
    rain_counter = rain_counter + 1
    IF (rain_counter == daily_timesteps) THEN
      !Re-set rain counter
      rain_counter = 0
      !Read next value of rainfall
      READ (090, *) rainfall
      READ (130, *) pet
    END IF

    !Terminate model run if steady-state has been reached
    !Note: 23.08.13. Add some on-screen confirmation here?
    IF (hydro_status == "steady") THEN
      EXIT
    END IF

    !On-screen counter to keep track of time
    WRITE(*, *) "t = ", time_counter

  END DO
  !End of hydrological loop

  !Write x and y indices to output files
  DO x = 1, x_extent
    DO y = 1, y_extent
      WRITE (110, '(1I10)') x
      WRITE (120, '(1I10)') y
    END DO
  END DO


END PROGRAM DigiBog_Hydro

!-------------------------------------------------------------------------------
! End of Program
!-------------------------------------------------------------------------------
