!-------------------------------------------------------------------------------
! cptecbam_crtm_mapping.f90
!-------------------------------------------------------------------------------
!BOI
! !TITLE: Surface Mapping Module for BAM Models
! !AUTHORS: João Gerd Zell de Mattos
! !AFFILIATION: CPTEC/INPE - Group on Data Assimilation Development (GDAD)
! !DATE: 2025-07-08
! !INTRODUCTION:
! This module provides the mapping interface between surface classification schemes
! used in CRTM (NPOESS/IGBP/USGS) and the vegetation/soil type definitions used in
! the BAM-SSiB and BAM-IBIS models. It defines default values for land, soil, and 
! vegetation types according to the CRTM User Guide (Table 4.11), and includes
! lookup tables for each of the supported classification systems.
!EOI


!-------------------------------------------------------------------------------
!BOP
! !MODULE: bam_surface_mapping.f90
! !DESCRIPTION:
!    This module provides mapping tables and routines to translate surface
!    classification schemes used by CRTM (NPOESS, IGBP, USGS, GFS) into the
!    corresponding codes expected by the land surface models BAM-SSiB and
!    BAM-IBIS. It defines:
!      1) Default CRTM surface type parameters.
!      2) Enumerations for each classification scheme (NPOESS, IGBP, USGS, GFS).
!      3) Bidirectional mapping arrays between SSiB/IBIS and CRTM schemes.
!      4) Zobler soil texture classification and mappings.
! Default values for CRTM_Surface_type are also defined to ensure robustness when mapping fails.
! These defaults follow the values suggested in Table 4.11 of the CRTM User Guide.

! !REVISION HISTORY:
!    07 Jul 2025 - J. G. de Mattos - Initial documentation.
!    08 Jul 2025 - J. G. Zell de Mattos - Initial version, fully documented using ProTeX.
! !SEE ALSO:
!  CRTM User Guide, Section 4.11.
!
! !INTERFACE:
module cptecbam_crtm_mapping
  use kinds,only: i_kind
  use crtm_module, only: crtm_irlandcoeff_classification, crtm_surface_type, microwave_sensor
  implicit none
  private
  
  external:: stop2

  public :: map_cptec_surface_types
  public :: map_bam_to_crtm
  public :: finalize_cptec_crtm_mapping

  !-----------------------------------------------------------------------!
  ! Final Maps to crtm
  !-----------------------------------------------------------------------!
  integer(i_kind),save, allocatable,dimension(:) :: map_to_crtm_ir
  integer(i_kind),save, allocatable,dimension(:) :: map_to_crtm_mw 
  integer(i_kind),save, allocatable,dimension(:) :: map_to_crtm_soil 

  !-----------------------------------------------------------------------!
  ! Default values for CRTM_Surface_type components (Table 4.11, CRTM UG) !
  !-----------------------------------------------------------------------!
  integer(i_kind), parameter :: DEFAULT_LAND_TYPE        = 1  ! compacted soil
  integer(i_kind), parameter :: DEFAULT_SOIL_TYPE        = 1  ! coarse loamy sand
  integer(i_kind), parameter :: DEFAULT_VEGETATION_TYPE  = 1  ! mixed/deciduous forest
  integer(i_kind), parameter :: DEFAULT_WATER_TYPE       = 1  ! sea water
  integer(i_kind), parameter :: DEFAULT_ICE_TYPE         = 1  ! new ice

  !----------------------------------------------------------------------------!
  ! Constants for land surface model identifiers                               !
  !----------------------------------------------------------------------------!
  integer(i_kind), parameter :: LSM_SSIB = 1
  integer(i_kind), parameter :: LSM_IBIS = 3

  !------------------------------------------------------------------------!
  ! Special flags for pure water and ice types                             !
  !------------------------------------------------------------------------!
  integer(i_kind), parameter :: WATER_TYPE = -1
  integer(i_kind), parameter :: ICE_TYPE   = -2


  !=============================================================================
  ! Esquema NPOESS (CRTM EmisCoeff) - Tipos de Superfície
  !=============================================================================
  ! Código NPOESS | Nome Simbólico                | Descrição
  !---------------+-------------------------------+------------------------------------------
  !             0 | INVALID_LAND                  | Inválido / sem dado
  !             1 | COMPACTED_SOIL                | Solo compactado
  !             2 | TILLED_SOIL                   | Solo arado
  !             3 | SAND                          | Areia
  !             4 | ROCK                          | Rocha
  !             5 | IRRIGATED_LOW_VEGETATION      | Vegetação baixa irrigada
  !             6 | MEADOW_GRASS                  | Gramíneas de campo
  !             7 | SCRUB                         | Arbustos
  !             8 | BROADLEAF_FOREST              | Floresta de folha larga
  !             9 | PINE_FOREST                   | Floresta de coníferas
  !            10 | TUNDRA                        | Tundra
  !            11 | GRASS_SOIL                    | Gramíneas com solo
  !            12 | BROADLEAF_PINE_FOREST         | Floresta mista (folha larga + pinheiro)
  !            13 | GRASS_SCRUB                   | Gramíneas e arbustos
  !            14 | SOIL_GRASS_SCRUB              | Solo, gramíneas e arbustos
  !            15 | URBAN_CONCRETE                | Área urbana/concreto
  !            16 | PINE_BRUSH                    | Pinheiros e arbustos
  !            17 | BROADLEAF_BRUSH               | Arbustos de folha larga
  !            18 | WET_SOIL                      | Solo úmido
  !            19 | SCRUB_SOIL                    | Arbustos e solo
  !            20 | BROADLEAF70_PINE30            | Mistura 70% folha larga + 30% pinheiro
  !=============================================================================
  integer(i_kind), parameter :: NPOESS_INVALID_LAND             =  0
  integer(i_kind), parameter :: NPOESS_COMPACTED_SOIL           =  1
  integer(i_kind), parameter :: NPOESS_TILLED_SOIL              =  2
  integer(i_kind), parameter :: NPOESS_SAND                     =  3
  integer(i_kind), parameter :: NPOESS_ROCK                     =  4
  integer(i_kind), parameter :: NPOESS_IRRIGATED_LOW_VEGETATION =  5
  integer(i_kind), parameter :: NPOESS_MEADOW_GRASS             =  6
  integer(i_kind), parameter :: NPOESS_SCRUB                    =  7
  integer(i_kind), parameter :: NPOESS_BROADLEAF_FOREST         =  8
  integer(i_kind), parameter :: NPOESS_PINE_FOREST              =  9
  integer(i_kind), parameter :: NPOESS_TUNDRA                   =  10
  integer(i_kind), parameter :: NPOESS_GRASS_SOIL               =  11
  integer(i_kind), parameter :: NPOESS_BROADLEAF_PINE_FOREST    =  12
  integer(i_kind), parameter :: NPOESS_GRASS_SCRUB              =  13
  integer(i_kind), parameter :: NPOESS_SOIL_GRASS_SCRUB         =  14
  integer(i_kind), parameter :: NPOESS_URBAN_CONCRETE           =  15
  integer(i_kind), parameter :: NPOESS_PINE_BRUSH               =  16
  integer(i_kind), parameter :: NPOESS_BROADLEAF_BRUSH          =  17
  integer(i_kind), parameter :: NPOESS_WET_SOIL                 =  18
  integer(i_kind), parameter :: NPOESS_SCRUB_SOIL               =  19
  integer(i_kind), parameter :: NPOESS_BROADLEAF70_PINE30       =  20

  !=============================================================================
  ! Esquema IGBP (CRTM EmisCoeff) - Tipos de Superfície
  !=============================================================================
  ! Código IGBP   | Nome Simbólico                | Descrição
  !---------------+-------------------------------+-----------------------------
  !             0 | NO_DATA                       | Sem dado / inválido
  !             1 | EVERGREEN_NEEDLELEAF_FOREST   | Floresta de coníferas perene
  !             2 | EVERGREEN_BROADLEAF_FOREST    | Floresta de folha larga perene
  !             3 | DECIDUOUS_NEEDLELEAF_FOREST   | Floresta de coníferas decídua (lariça)
  !             4 | DECIDUOUS_BROADLEAF_FOREST    | Floresta de folha larga decídua
  !             5 | MIXED_FOREST                  | Florestas mistas (needle + broadleaf)
  !             6 | CLOSED_SHRUBLANDS             | Arbustais fechados
  !             7 | OPEN_SHRUBLANDS               | Arbustais abertos
  !             8 | WOODY_SAVANNAS                | Savanas arbóreas
  !             9 | SAVANNAS                      | Savanas gramíneas
  !            10 | GRASSLANDS                    | Campos / pradarias
  !            11 | PERMANENT_WETLANDS            | Áreas alagadas permanentes
  !            12 | CROPLANDS                     | Áreas agrícolas
  !            13 | URBAN_AND_BUILT_UP            | Áreas urbanas e edificadas
  !            14 | CROPLAND_NATURAL_VEG_MOSAIC   | Mosaico agrícola e vegetação natural
  !            15 | SNOW_AND_ICE                  | Neve e gelo
  !            16 | BARREN_OR_SPARSE_VEGETATION   | Solo exposto ou vegetação esparsa
  !            17 | WATER_BODIES                  | Corpos d’água
  !=============================================================================
  integer(i_kind), parameter :: IGBP_NO_DATA                     =  0
  integer(i_kind), parameter :: IGBP_EVERGREEN_NEEDLELEAF_FOREST =  1
  integer(i_kind), parameter :: IGBP_EVERGREEN_BROADLEAF_FOREST  =  2
  integer(i_kind), parameter :: IGBP_DECIDUOUS_NEEDLELEAF_FOREST =  3
  integer(i_kind), parameter :: IGBP_DECIDUOUS_BROADLEAF_FOREST  =  4
  integer(i_kind), parameter :: IGBP_MIXED_FOREST                =  5
  integer(i_kind), parameter :: IGBP_CLOSED_SHRUBLANDS           =  6
  integer(i_kind), parameter :: IGBP_OPEN_SHRUBLANDS             =  7
  integer(i_kind), parameter :: IGBP_WOODY_SAVANNAS              =  8
  integer(i_kind), parameter :: IGBP_SAVANNAS                    =  9
  integer(i_kind), parameter :: IGBP_GRASSLANDS                  = 10
  integer(i_kind), parameter :: IGBP_PERMANENT_WETLANDS          = 11
  integer(i_kind), parameter :: IGBP_CROPLANDS                   = 12
  integer(i_kind), parameter :: IGBP_URBAN_AND_BUILT_UP          = 13
  integer(i_kind), parameter :: IGBP_CROPLAND_NATURAL_VEG_MOSAIC = 14
  integer(i_kind), parameter :: IGBP_SNOW_AND_ICE                = 15
  integer(i_kind), parameter :: IGBP_BARREN_OR_SPARSE_VEGETATION = 16
  integer(i_kind), parameter :: IGBP_WATER_BODIES                = 17

  !================================================================================
  ! Esquema USGS (CRTM EmisCoeff) - Tipos de Superfície
  !================================================================================
  ! Código USGS   | Nome Simbólico                   | Descrição
  !---------------+----------------------------------+-----------------------------
  !             1 | URBAN_BUILTUP_LAND               | urban and built-up land
  !             2 | DRYLAND_CROPLAND_PASTURE         | dryland cropland and pasture
  !             3 | IRRIGATED_CROPLAND_PASTURE       | irrigated cropland and pasture
  !             4 | MIXED_DRY_IRRIG_CROPLAND_PASTURE | mixed dryland/irrigated cropland and pasture
  !             5 | CROPLAND_GRASSLAND_MOSAIC        | cropland/grassland mosaic
  !             6 | CROPLAND_WOODLAND_MOSAIC         | cropland/woodland mosaic
  !             7 | GRASSLAND                        | grassland
  !             8 | SHRUBLAND                        | shrubland
  !             9 | MIXED_SHRUB_GRASSLAND            | mixed shrubland/grassland
  !            10 | SAVANNA                          | savanna
  !            11 | DECIDUOUS_BROADLEAF_FOREST       | deciduous broadleaf forest
  !            12 | DECIDUOUS_NEEDLELEAF_FOREST      | deciduous needleleaf forest
  !            13 | EVERGREEN_BROADLEAF_FOREST       | evergreen broadleaf forest
  !            14 | EVERGREEN_NEEDLELEAF_FOREST      | evergreen needleleaf forest
  !            15 | MIXED_FOREST                     | mixed forest
  !            16 | WATER_BODIES                     | water bodies
  !            17 | HERBACEOUS_WETLAND               | herbaceous wetland
  !            18 | WOODED_WETLAND                   | wooded wetland
  !            19 | BARREN_OR_SPARSE_VEGETATED       | barren or sparsely vegetated
  !            20 | HERBACEOUS_TUNDRA                | herbaceous tundra
  !            21 | WOODED_TUNDRA                    | wooded tundra
  !            22 | MIXED_TUNDRA                     | mixed tundra
  !            23 | BARE_GROUND_TUNDRA               | bare ground tundra
  !            24 | SNOW_OR_ICE                      | snow or ice
  !            25 | PLAYA                            | playa
  !            26 | LAVA                             | lava
  !            27 | WHITE_SAND                       | white sand
  !================================================================================
  integer(i_kind), parameter :: USGS_URBAN_BUILTUP_LAND               =  1
  integer(i_kind), parameter :: USGS_DRYLAND_CROPLAND_PASTURE         =  2
  integer(i_kind), parameter :: USGS_IRRIGATED_CROPLAND_PASTURE       =  3
  integer(i_kind), parameter :: USGS_MIXED_DRY_IRRIG_CROPLAND_PASTURE =  4
  integer(i_kind), parameter :: USGS_CROPLAND_GRASSLAND_MOSAIC        =  5
  integer(i_kind), parameter :: USGS_CROPLAND_WOODLAND_MOSAIC         =  6
  integer(i_kind), parameter :: USGS_GRASSLAND                        =  7
  integer(i_kind), parameter :: USGS_SHRUBLAND                        =  8
  integer(i_kind), parameter :: USGS_MIXED_SHRUB_GRASSLAND            =  9
  integer(i_kind), parameter :: USGS_SAVANNA                          = 10
  integer(i_kind), parameter :: USGS_DECIDUOUS_BROADLEAF_FOREST       = 11
  integer(i_kind), parameter :: USGS_DECIDUOUS_NEEDLELEAF_FOREST      = 12
  integer(i_kind), parameter :: USGS_EVERGREEN_BROADLEAF_FOREST       = 13
  integer(i_kind), parameter :: USGS_EVERGREEN_NEEDLELEAF_FOREST      = 14
  integer(i_kind), parameter :: USGS_MIXED_FOREST                     = 15
  integer(i_kind), parameter :: USGS_WATER_BODIES                     = 16
  integer(i_kind), parameter :: USGS_HERBACEOUS_WETLAND               = 17
  integer(i_kind), parameter :: USGS_WOODED_WETLAND                   = 18
  integer(i_kind), parameter :: USGS_BARREN_OR_SPARSE_VEGETATED       = 19
  integer(i_kind), parameter :: USGS_HERBACEOUS_TUNDRA                = 20
  integer(i_kind), parameter :: USGS_WOODED_TUNDRA                    = 21
  integer(i_kind), parameter :: USGS_MIXED_TUNDRA                     = 22
  integer(i_kind), parameter :: USGS_BARE_GROUND_TUNDRA               = 23
  integer(i_kind), parameter :: USGS_SNOW_OR_ICE                      = 24
  integer(i_kind), parameter :: USGS_PLAYA                            = 25
  integer(i_kind), parameter :: USGS_LAVA                             = 26
  integer(i_kind), parameter :: USGS_WHITE_SAND                       = 27

  !=============================================================================
  ! Esquema GFS (CRTM EmisCoeff) - Tipos de Superfície
  !=============================================================================
  ! Código GFS    | Nome Simbólico                | Descrição
  !---------------+-------------------------------+-----------------------------
  !             0 | WATER                         | water
  !             1 | BROADLEAF_EVERGREEN           | broadleaf-evergreen (tropical forest)
  !             2 | BROAD_DECIDUOUS               | broad-deciduous trees
  !             3 | MIXED_FOREST                  | broadleaf & needleleaf (mixed forest)
  !             4 | NEEDLELEAF_EVERGREEN          | needleleaf-evergreen trees
  !             5 | NEEDLELEAF_DECIDUOUS          | needleleaf-deciduous trees (larch)
  !             6 | SAVANNA                       | broadleaf trees with ground cover (savanna)
  !             7 | GRASSLAND                     | ground cover only (perennial)
  !             8 | SHRUB_COVER                   | broad leaf shrubs w/ ground cover
  !             9 | SHRUB_BARE                    | broadleaf shrubs with bare soil
  !            10 | TUNDRA_SHRUBS                 | dwarf trees & shrubs w/ ground cover (tundra)
  !            11 | BARE_SOIL                     | bare soil
  !            12 | CROPLANDS                     | crops
  !            13 | PERMANENT_ICE                 | permanent ice
  !=============================================================================
  integer(i_kind), parameter :: GFS_WATER                 =  0  ! water
  integer(i_kind), parameter :: GFS_BROADLEAF_EVERGREEN   =  1  ! broadleaf-evergreen (tropical forest)
  integer(i_kind), parameter :: GFS_BROAD_DECIDUOUS       =  2  ! broad-deciduous trees
  integer(i_kind), parameter :: GFS_MIXED_FOREST          =  3  ! broadleaf & needleleaf (mixed forest)
  integer(i_kind), parameter :: GFS_NEEDLELEAF_EVERGREEN  =  4  ! needleleaf-evergreen trees
  integer(i_kind), parameter :: GFS_NEEDLELEAF_DECIDUOUS  =  5  ! needleleaf-deciduous trees (larch)
  integer(i_kind), parameter :: GFS_SAVANNA               =  6  ! broadleaf trees with ground cover (savanna)
  integer(i_kind), parameter :: GFS_GRASSLAND             =  7  ! ground cover only (perennial)
  integer(i_kind), parameter :: GFS_SHRUB_COVER           =  8  ! broad leaf shrubs w/ ground cover
  integer(i_kind), parameter :: GFS_SHRUB_BARE            =  9  ! broadleaf shrubs with bare soil
  integer(i_kind), parameter :: GFS_TUNDRA_SHRUBS         = 10  ! dwarf trees & shrubs w/ ground cover (tundra)
  integer(i_kind), parameter :: GFS_BARE_SOIL             = 11  ! bare soil
  integer(i_kind), parameter :: GFS_CROPLANDS             = 12  ! crops
  integer(i_kind), parameter :: GFS_PERMANENT_ICE         = 13  ! permanent ice

  !=============================================================================
  ! Number of vegetation types in each LSM
  !=============================================================================
  integer(i_kind), parameter :: SSIB_VEG_N = 13    ! SSiB 0..13
  integer(i_kind), parameter :: IBIS_VEG_N = 15    ! IBIS 0..15

  !=============================================================================
  ! Tabela de equivalência SSiB → CRTM NPOESS (Infrared/VIS)
  !
  ! Propósito:
  !   Converter as classes de vegetação do modelo SSiB para os códigos de
  !   superfície do esquema NPOESS usados no CRTM.
  !
  ! Critérios de mapeamento:
  ! - Broadleaf-evergreen (SSiB 1) → BROADLEAF_FOREST
  ! - Broadleaf-deciduous e Mixed forest (SSiB 2,3) → BROADLEAF_PINE_FOREST
  ! - Needleleaf-evergreen e Deciduous needleleaf (SSiB 4,5) → PINE_FOREST
  ! - Savana (SSiB 6) → GRASS_SCRUB (mistura gramíneas + arbustos)
  ! - Grassland (SSiB 7) → MEADOW_GRASS
  ! - Shrublands (SSiB 8,9) → SCRUB / SCRUB_SOIL conforme solo exposto
  ! - Tundra (SSiB 10) → TUNDRA
  ! - Barren soil (SSiB 11) → COMPACTED_SOIL
  ! - Crops (SSiB 12) → IRRIGATED_LOW_VEGETATION
  ! - Água (SSiB 0) → WATER_TYPE; Gelo (SSiB 13) → ICE_TYPE
  !=============================================================================
  integer(i_kind), dimension(0:SSIB_VEG_N), parameter :: ssib_to_npoess = (/ &
     WATER_TYPE,                     & ! 0: water
     NPOESS_BROADLEAF_FOREST,        & ! 1: broadleaf-evergreen
     NPOESS_BROADLEAF_PINE_FOREST,   & ! 2: broad-deciduous
     NPOESS_BROADLEAF_PINE_FOREST,   & ! 3: mixed forest
     NPOESS_PINE_FOREST,             & ! 4: needleleaf-evergreen
     NPOESS_PINE_FOREST,             & ! 5: needleleaf‐deciduous trees (larch)
     NPOESS_GRASS_SCRUB,             & ! 6: broadleaf trees with ground cover (savanna)
     NPOESS_MEADOW_GRASS,            & ! 7: ground cover only (perennial)
     NPOESS_GRASS_SCRUB,             & ! 8: broad leaf shrubs w/ ground cover
     NPOESS_SCRUB_SOIL,              & ! 9: broadleaf shrubs with bare soil
     NPOESS_TUNDRA,                  & !10: dwarf trees & shrubs w/ground cover (tundra)
     NPOESS_COMPACTED_SOIL,          & !11: bare soil
     NPOESS_IRRIGATED_LOW_VEGETATION,& !12: crops
     ICE_TYPE                        & !13: ice
  /)

  !=============================================================================
  ! Tabela de equivalência IBIS → CRTM NPOESS (Infrared/VIS)
  !
  ! Propósito:
  !   Mapear as classes de vegetação do modelo IBIS para os códigos de
  !   superfície do esquema NPOESS usados no CRTM.
  !
  ! Critérios de mapeamento:
  ! - Tropical evergreen e deciduous (IBIS 1,2), Temperate broadleaf (IBIS 3,5) →
  !     BROADLEAF_FOREST
  ! - Evergreen conifer e Boreal (IBIS 4,6,7) → PINE_FOREST
  ! - Mixed woodland/forest (IBIS 8) → BROADLEAF_PINE_FOREST
  ! - Savana e Grassland/steppe (IBIS 9,10) → MEADOW_GRASS
  ! - Shrublands dense/open (IBIS 11,12) → SCRUB / GRASS_SCRUB
  ! - Tundra (IBIS 13) → TUNDRA
  ! - Desert (IBIS 14) → COMPACTED_SOIL
  ! - Água (IBIS 0) → WATER_TYPE; Polar desert/ice (IBIS 15) → ICE_TYPE
  !=============================================================================
  integer(i_kind), dimension(0:IBIS_VEG_N), parameter :: ibis_to_npoess = (/ &
     WATER_TYPE,                     & !  0: Ocean/Lakes/Rivers       -> (via Water_Type)
     NPOESS_BROADLEAF_FOREST,        & !  1: tropical evergreen
     NPOESS_BROADLEAF_PINE_FOREST,   & !  2: tropical deciduous
     NPOESS_BROADLEAF_FOREST,        & !  3: temperate evergreen broadleaf
     NPOESS_PINE_FOREST,             & !  4: temperate evergreen conifer
     NPOESS_BROADLEAF_PINE_FOREST,   & !  5: temperate deciduous forest
     NPOESS_PINE_FOREST,             & !  6: boreal evergreen
     NPOESS_PINE_FOREST,             & !  7: boreal deciduous
     NPOESS_BROADLEAF_PINE_FOREST,   & !  8: mixed forest/woodland
     NPOESS_MEADOW_GRASS,            & !  9: Savana
     NPOESS_MEADOW_GRASS,            & ! 10: grassland/steppe
     NPOESS_SCRUB,                   & ! 11: dense shrubland
     NPOESS_GRASS_SCRUB,             & ! 12: open shrubland
     NPOESS_TUNDRA,                  & ! 13: Tundra
     NPOESS_COMPACTED_SOIL,          & ! 14: Desert
     ICE_TYPE                        & ! 15: polar desert/rock/ice    -> (via Ice_Type)
  /)

  !=============================================================================
  ! Tabela de equivalência SSiB → CRTM IGBP (Infrared/VIS)
  !
  ! Propósito:
  !   Traduzir as classes de vegetação do modelo SSiB para os códigos IGBP
  !   usados no CRTM.
  !
  ! Critérios de mapeamento:
  ! - Broadleaf-evergreen (SSiB 1) → EVERGREEN_BROADLEAF_FOREST
  ! - Broadleaf-deciduous (SSiB 2) → DECIDUOUS_BROADLEAF_FOREST
  ! - Mixed forest (SSiB 3) → MIXED_FOREST
  ! - Needleleaf-evergreen (SSiB 4) → EVERGREEN_NEEDLELEAF_FOREST
  ! - Deciduous needleleaf (SSiB 5) → DECIDUOUS_NEEDLELEAF_FOREST
  ! - Savana (SSiB 6) → SAVANNAS
  ! - Grassland (SSiB 7) → GRASSLANDS
  ! - Shrublands (SSiB 8,9) → CLOSED_SHRUBLANDS / OPEN_SHRUBLANDS
  ! - Tundra (SSiB 10) → BARREN_OR_SPARSE_VEGETATION (tundra genérica)
  ! - Bare soil (SSiB 11) → BARREN_OR_SPARSE_VEGETATION
  ! - Crops (SSiB 12) → CROPLANDS
  ! - Água (SSiB 0) → WATER_BODIES; Gelo (SSiB 13) → SNOW_AND_ICE
  !=============================================================================
  integer(i_kind), dimension(0:SSIB_VEG_N), parameter :: ssib_to_igbp = (/ &
     WATER_TYPE,                             & !  0: water
     IGBP_EVERGREEN_BROADLEAF_FOREST,        & !  1: broadleaf-evergreen
     IGBP_DECIDUOUS_BROADLEAF_FOREST,        & !  2: broad-deciduous
     IGBP_MIXED_FOREST,                      & !  3: mixed forest
     IGBP_EVERGREEN_NEEDLELEAF_FOREST,       & !  4: needleleaf-evergreen
     IGBP_DECIDUOUS_NEEDLELEAF_FOREST,       & !  5: needleleaf-deciduous
     IGBP_SAVANNAS,                          & !  6: savanna
     IGBP_GRASSLANDS,                        & !  7: ground cover only
     IGBP_CLOSED_SHRUBLANDS,                 & !  8: woody shrubs
     IGBP_OPEN_SHRUBLANDS,                   & !  9: shrubs with bare soil
     IGBP_BARREN_OR_SPARSE_VEGETATION,       & ! 10: tundra → barren
     IGBP_BARREN_OR_SPARSE_VEGETATION,       & ! 11: bare soil
     IGBP_CROPLANDS,                         & ! 12: crops
     ICE_TYPE                                & ! 13: ice
  /)

  !=============================================================================
  ! Tabela de equivalência IBIS → CRTM IGBP (Infrared/VIS)
  !
  ! Propósito:
  !   Converter classes IBIS em códigos IGBP para uso no CRTM.
  !
  ! Critérios de mapeamento:
  ! - Tropical & Temperate broadleaf (IBIS 1,2,3,5) → EVERGREEN/DECIDUOUS_BROADLEAF_FOREST
  ! - Coniferous evergreen/deciduous (IBIS 4,6,7) → EVERGREEN/DECIDUOUS_NEEDLELEAF_FOREST
  ! - Mixed forest/woodland (IBIS 8) → MIXED_FOREST
  ! - Savana (IBIS 9) → SAVANNAS
  ! - Grassland/steppe (IBIS 10) → GRASSLANDS
  ! - Shrublands (IBIS 11,12) → CLOSED_SHRUBLANDS / OPEN_SHRUBLANDS
  ! - Tundra & Desert (IBIS 13,14) → BARREN_OR_SPARSE_VEGETATION
  ! - Água (IBIS 0) → WATER_BODIES; Polar desert/ice (IBIS 15) → SNOW_AND_ICE
  !=============================================================================
  integer(i_kind), dimension(0:IBIS_VEG_N), parameter :: ibis_to_igbp = (/ &
     WATER_TYPE,                             & !  0: water
     IGBP_EVERGREEN_BROADLEAF_FOREST,        & !  1: tropical evergreen
     IGBP_DECIDUOUS_BROADLEAF_FOREST,        & !  2: tropical deciduous
     IGBP_EVERGREEN_BROADLEAF_FOREST,        & !  3: temperate evergreen broadleaf
     IGBP_EVERGREEN_NEEDLELEAF_FOREST,       & !  4: temperate evergreen conifer
     IGBP_MIXED_FOREST,                      & !  5: temperate deciduous forest
     IGBP_EVERGREEN_NEEDLELEAF_FOREST,       & !  6: boreal evergreen
     IGBP_DECIDUOUS_NEEDLELEAF_FOREST,       & !  7: boreal deciduous
     IGBP_MIXED_FOREST,                      & !  8: mixed forest/woodland
     IGBP_SAVANNAS,                          & !  9: savanna
     IGBP_GRASSLANDS,                        & ! 10: grassland/steppe
     IGBP_CLOSED_SHRUBLANDS,                 & ! 11: dense shrubland
     IGBP_OPEN_SHRUBLANDS,                   & ! 12: open shrubland
     IGBP_BARREN_OR_SPARSE_VEGETATION,       & ! 13: tundra → barren
     IGBP_BARREN_OR_SPARSE_VEGETATION,       & ! 14: desert → barren
     ICE_TYPE                                & ! 15: ice
  /)

  !=============================================================================
  ! Tabela de equivalência SSiB → CRTM USGS (Infrared/VIS)
  !
  ! Propósito:
  !   Mapear classes de vegetação SSiB para códigos USGS no CRTM.
  !
  ! Critérios de mapeamento:
  ! - Broadleaf-evergreen/deciduous (SSiB 1,2) → EVERGREEN/DECIDUOUS_BROADLEAF_FOREST
  ! - Mixed forest (SSiB 3) → MIXED_FOREST
  ! - Needleleaf (SSiB 4,5) → EVERGREEN/DECIDUOUS_NEEDLELEAF_FOREST
  ! - Savana (SSiB 6) → SAVANNA
  ! - Grassland (SSiB 7) → GRASSLAND
  ! - Shrublands (SSiB 8,9) → SHRUBLAND / MIXED_SHRUB_GRASSLAND
  ! - Tundra (SSiB 10) → MIXED_TUNDRA
  ! - Bare soil (SSiB 11) → BARREN_OR_SPARSE_VEGETATED
  ! - Crops (SSiB 12) → DRYLAND_CROPLAND_PASTURE
  ! - Água (SSiB 0) → WATER_BODIES; Gelo (SSiB 13) → SNOW_OR_ICE
  !=============================================================================
  integer(i_kind), dimension(0:SSIB_VEG_N), parameter :: ssib_to_usgs = (/ &
     USGS_WATER_BODIES,                      & !  0: water
     USGS_EVERGREEN_BROADLEAF_FOREST,        & !  1: broadleaf-evergreen
     USGS_DECIDUOUS_BROADLEAF_FOREST,        & !  2: broad-deciduous
     USGS_MIXED_FOREST,                      & !  3: mixed forest
     USGS_EVERGREEN_NEEDLELEAF_FOREST,       & !  4: needleleaf-evergreen
     USGS_DECIDUOUS_NEEDLELEAF_FOREST,       & !  5: needleleaf-deciduous
     USGS_SAVANNA,                           & !  6: savanna
     USGS_GRASSLAND,                         & !  7: ground cover only
     USGS_SHRUBLAND,                         & !  8: woody shrubs
     USGS_MIXED_SHRUB_GRASSLAND,             & !  9: shrubs with soil
     USGS_MIXED_TUNDRA,                      & ! 10: tundra
     USGS_BARREN_OR_SPARSE_VEGETATED,        & ! 11: bare soil
     USGS_DRYLAND_CROPLAND_PASTURE,          & ! 12: crops
     USGS_SNOW_OR_ICE                        & ! 13: ice
  /)

  !=============================================================================
  ! Tabela de equivalência IBIS → CRTM USGS (Infrared/VIS)
  !
  ! Propósito:
  !   Converter classes IBIS em códigos USGS para o CRTM.
  !
  ! Critérios de mapeamento:
  ! - Broadleaf forests (IBIS 1,2,3,5) → EVERGREEN/DECIDUOUS_BROADLEAF_FOREST
  ! - Conifer forests (IBIS 4,6,7) → EVERGREEN/DECIDUOUS_NEEDLELEAF_FOREST
  ! - Mixed forest/woodland (IBIS 8) → MIXED_FOREST
  ! - Savana & Grassland (IBIS 9,10) → SAVANNA / GRASSLAND
  ! - Shrublands (IBIS 11,12) → SHRUBLAND / MIXED_SHRUB_GRASSLAND
  ! - Tundra (IBIS 13) → MIXED_TUNDRA; Desert (IBIS 14) → BARREN_OR_SPARSE_VEGETATED
  ! - Água (IBIS 0) → WATER_BODIES; Polar desert/ice (IBIS 15) → SNOW_OR_ICE
  !=============================================================================
  integer(i_kind), dimension(0:IBIS_VEG_N), parameter :: ibis_to_usgs = (/ &
     USGS_WATER_BODIES,                      & !  0: water
     USGS_EVERGREEN_BROADLEAF_FOREST,        & !  1: tropical evergreen
     USGS_DECIDUOUS_BROADLEAF_FOREST,        & !  2: tropical deciduous
     USGS_EVERGREEN_BROADLEAF_FOREST,        & !  3: temperate evergreen broadleaf
     USGS_EVERGREEN_NEEDLELEAF_FOREST,       & !  4: temperate evergreen conifer
     USGS_MIXED_FOREST,                      & !  5: temperate deciduous forest
     USGS_EVERGREEN_NEEDLELEAF_FOREST,       & !  6: boreal evergreen
     USGS_DECIDUOUS_NEEDLELEAF_FOREST,       & !  7: boreal deciduous
     USGS_MIXED_FOREST,                      & !  8: mixed forest/woodland
     USGS_SAVANNA,                           & !  9: savanna
     USGS_GRASSLAND,                         & ! 10: grassland/steppe
     USGS_SHRUBLAND,                         & ! 11: dense shrubland
     USGS_SHRUBLAND,                         & ! 12: open shrubland
     USGS_MIXED_TUNDRA,                      & ! 13: tundra
     USGS_BARREN_OR_SPARSE_VEGETATED,        & ! 14: desert → barren
     USGS_SNOW_OR_ICE                        & ! 15: ice
  /)

  !=============================================================================
  ! Tabela de equivalência SSiB → GFS
  !
  ! Propósito:
  !   Mapear as classes de vegetação do modelo SSiB para os códigos GFS.
  !
  ! Critérios de mapeamento:
  ! - Florestas broadleaf-evergreen (1) → GFS_BROADLEAF_EVERGREEN
  ! - Florestas broadleaf-deciduous (2) → GFS_BROAD_DECIDUOUS
  ! - Florestas mistas (3) → GFS_MIXED_FOREST
  ! - Florestas coníferas perenes (4) → GFS_NEEDLELEAF_EVERGREEN
  ! - Lariças (5) → GFS_NEEDLELEAF_DECIDUOUS
  ! - Savana (6) → GFS_SAVANNA
  ! - Gramíneas puras (7) → GFS_GRASSLAND
  ! - Arbustos com cobertura (8) → GFS_SHRUB_COVER
  ! - Arbustos com solo exposto (9) → GFS_SHRUB_BARE
  ! - Tundra (10) → GFS_TUNDRA_SHRUBS
  ! - Solo exposto (11) → GFS_BARE_SOIL
  ! - Cultivos (12) → GFS_CROPLANDS
  ! - Água (0) → WATER_TYPE; Gelo (13) → ICE_TYPE
  !=============================================================================
  integer(i_kind), dimension(0:SSIB_VEG_N), parameter :: ssib_to_gfs = (/ &
     WATER_TYPE,                   & !  0: water
     GFS_BROADLEAF_EVERGREEN,      & !  1: broadleaf-evergreen
     GFS_BROAD_DECIDUOUS,          & !  2: broad-deciduous
     GFS_MIXED_FOREST,             & !  3: mixed forest
     GFS_NEEDLELEAF_EVERGREEN,     & !  4: needleleaf-evergreen
     GFS_NEEDLELEAF_DECIDUOUS,     & !  5: needleleaf-deciduous
     GFS_SAVANNA,                  & !  6: savanna
     GFS_GRASSLAND,                & !  7: ground cover only
     GFS_SHRUB_COVER,              & !  8: shrubs w/ ground cover
     GFS_SHRUB_BARE,               & !  9: shrubs with bare soil
     GFS_TUNDRA_SHRUBS,            & ! 10: dwarf shrubs (tundra)
     GFS_BARE_SOIL,                & ! 11: bare soil
     GFS_CROPLANDS,                & ! 12: crops
     ICE_TYPE                      & ! 13: permanent ice
  /)

  !=============================================================================
  ! Tabela de equivalência IBIS → GFS
  !
  ! Propósito:
  !   Mapear as classes de vegetação do modelo IBIS para os códigos GFS.
  !
  ! Critérios de mapeamento:
  ! - Tropical evergreen (1) → GFS_BROADLEAF_EVERGREEN
  ! - Tropical deciduous (2) → GFS_BROAD_DECIDUOUS
  ! - Temperate mixed & woodland (3,8) → GFS_MIXED_FOREST
  ! - Coníferas perenes e boreais (4,6) → GFS_NEEDLELEAF_EVERGREEN
  ! - Lariças (7) → GFS_NEEDLELEAF_DECIDUOUS
  ! - Savana (9) → GFS_SAVANNA
  ! - Grassland/steppe (10) → GFS_GRASSLAND
  ! - Dense/open shrubland (11,12) → GFS_SHRUB_COVER / GFS_SHRUB_BARE
  ! - Tundra (13) → GFS_TUNDRA_SHRUBS
  ! - Deserto (14) → GFS_BARE_SOIL
  ! - Água (0) → WATER_TYPE; Polar desert/ice (15) → ICE_TYPE
  !=============================================================================
  integer(i_kind), dimension(0:IBIS_VEG_N), parameter :: ibis_to_gfs = (/ &
     WATER_TYPE,                   & !  0: water
     GFS_BROADLEAF_EVERGREEN,      & !  1: tropical evergreen
     GFS_BROAD_DECIDUOUS,          & !  2: tropical deciduous
     GFS_MIXED_FOREST,             & !  3: mixed forest
     GFS_NEEDLELEAF_EVERGREEN,     & !  4: temperate evergreen conifer
     GFS_MIXED_FOREST,             & !  5: temperate deciduous forest
     GFS_NEEDLELEAF_EVERGREEN,     & !  6: boreal evergreen
     GFS_NEEDLELEAF_DECIDUOUS,     & !  7: boreal deciduous
     GFS_MIXED_FOREST,             & !  8: mixed forest/woodland
     GFS_SAVANNA,                  & !  9: savanna
     GFS_GRASSLAND,                & ! 10: grassland/steppe
     GFS_SHRUB_COVER,              & ! 11: dense shrubland
     GFS_SHRUB_BARE,               & ! 12: open shrubland
     GFS_TUNDRA_SHRUBS,            & ! 13: tundra
     GFS_BARE_SOIL,                & ! 14: desert → bare soil
     ICE_TYPE                      & ! 15: ice
  /)

  !=============================================================================
  ! Flags especiais
  !=============================================================================
  integer(i_kind), parameter :: LAND_ICE    =  9   ! glacial land ice
  integer(i_kind), parameter :: FARMLAND    =  8   ! organic / farmland

  !=============================================================================
  ! Zobler Soil Type Classification Scheme
  !=============================================================================
  ! Index : Texture         : Descrição
  !-----------------------------------------------------!
  !  1   : coarse           : loamy sand
  !  2   : medium           : silty clay loam
  !  3   : fine             : light clay
  !  4   : coarse-medium    : sandy loam
  !  5   : coarse-fine      : sandy clay
  !  6   : medium-fine      : clay loam
  !  7   : coarse-med-fine  : sand clay loam
  !  8   : organic          : farmland
  !  9   : glacial land ice : ice over land
  !=============================================================================
  integer(i_kind), parameter :: ZOBLER_COARSE            = 1   ! loamy sand
  integer(i_kind), parameter :: ZOBLER_MEDIUM            = 2   ! silty clay loam
  integer(i_kind), parameter :: ZOBLER_FINE              = 3   ! light clay
  integer(i_kind), parameter :: ZOBLER_COARSE_MEDIUM     = 4   ! sandy loam
  integer(i_kind), parameter :: ZOBLER_COARSE_FINE       = 5   ! sandy clay
  integer(i_kind), parameter :: ZOBLER_MEDIUM_FINE       = 6   ! clay loam
  integer(i_kind), parameter :: ZOBLER_COARSE_MED_FINE   = 7   ! sand clay loam
  integer(i_kind), parameter :: ZOBLER_ORGANIC           = 8   ! farmland (organic)
  integer(i_kind), parameter :: ZOBLER_GLACIAL           = 9   ! ice over land

 !=============================================================================
  ! Número de tipos de solo em cada LSM
  !=============================================================================
  integer(i_kind), parameter :: SSIB_SOIL_N = 12   ! SSiB classes 1..12
  integer(i_kind), parameter :: IBIS_SOIL_N = 12   ! IBIS classes 1..12

  !=============================================================================
  ! Tabela de equivalência SSiB → Zobler
  !
  ! Propósito:
  !   Converter as classes de solo do SSiB para a classificação Zobler.
  !
  ! Critérios de mapeamento:
  ! - Sand (1) → ZOBLER_COARSE_MEDIUM (sandy loam), pois SSiB pure sand não existe em Zobler
  ! - Loamy Sand (2) → ZOBLER_COARSE (loamy sand)
  ! - Sandy Loam (3) → ZOBLER_COARSE_MEDIUM (sandy loam)
  ! - Silt Loam (4) → ZOBLER_MEDIUM_FINE (clay loam), aproximação por textura fina
  ! - Loam (5) → ZOBLER_MEDIUM_FINE (clay loam), textura intermediária
  ! - Sandy Clay Loam (6) → ZOBLER_COARSE_MED_FINE (sand clay loam)
  ! - Silty Clay Loam (7) → ZOBLER_MEDIUM (silty clay loam)
  ! - Clay Loam (8) → ZOBLER_MEDIUM_FINE (clay loam)
  ! - Sandy Clay (9) → ZOBLER_COARSE_FINE (sandy clay)
  ! - Silty Clay (10) → ZOBLER_FINE (light clay), aproximação para solo muito fino
  ! - Clay (11) → ZOBLER_FINE (light clay)
  ! - Silt (12) → ZOBLER_MEDIUM (silty clay loam), aproximação por textura média
  !=============================================================================
  integer(i_kind), dimension(1:SSIB_SOIL_N), parameter :: ssib_to_zobler = (/ &
     ZOBLER_COARSE_MEDIUM,  & ! 1: Sand
     ZOBLER_COARSE,         & ! 2: Loamy Sand
     ZOBLER_COARSE_MEDIUM,  & ! 3: Sandy Loam
     ZOBLER_MEDIUM_FINE,    & ! 4: Silt Loam
     ZOBLER_MEDIUM_FINE,    & ! 5: Loam
     ZOBLER_COARSE_MED_FINE,& ! 6: Sandy Clay Loam
     ZOBLER_MEDIUM,         & ! 7: Silty Clay Loam
     ZOBLER_MEDIUM_FINE,    & ! 8: Clay Loam
     ZOBLER_COARSE_FINE,    & ! 9: Sandy Clay
     ZOBLER_FINE,           & !10: Silty Clay
     ZOBLER_FINE,           & !11: Clay
     ZOBLER_MEDIUM          & !12: Silt
  /)

  !=============================================================================
  ! Tabela de equivalência IBIS → Zobler
  !
  ! Propósito:
  !   Converter as classes de solo do IBIS para a classificação Zobler.
  !
  ! Critérios de mapeamento:
  ! - Sand (1) → ZOBLER_COARSE_MEDIUM (sandy loam)
  ! - Loamy Sand (2) → ZOBLER_COARSE (loamy sand)
  ! - Sand Loam (3) → ZOBLER_COARSE_MEDIUM (sandy loam)
  ! - Loam (4) → ZOBLER_MEDIUM_FINE (clay loam)
  ! - Silty Loam (5) → ZOBLER_MEDIUM_FINE (clay loam)
  ! - Sand Clay Loam (6) → ZOBLER_COARSE_MED_FINE (sand clay loam)
  ! - Clay Loam (7) → ZOBLER_MEDIUM_FINE (clay loam)
  ! - Silty Clay Loam (8) → ZOBLER_MEDIUM (silty clay loam)
  ! - Sandy Clay (9) → ZOBLER_COARSE_FINE (sandy clay)
  ! - Silty Clay (10) → ZOBLER_FINE (light clay)
  ! - Clay (11) → ZOBLER_FINE (light clay)
  ! - Organic (12) → ZOBLER_ORGANIC (farmland)
  !=============================================================================
  integer(i_kind), dimension(1:IBIS_SOIL_N), parameter :: ibis_to_zobler = (/ &
     ZOBLER_COARSE_MEDIUM,  & ! 1: Sand
     ZOBLER_COARSE,         & ! 2: Loamy Sand
     ZOBLER_COARSE_MEDIUM,  & ! 3: Sand Loam
     ZOBLER_MEDIUM_FINE,    & ! 4: Loam
     ZOBLER_MEDIUM_FINE,    & ! 5: Silty Loam
     ZOBLER_COARSE_MED_FINE,& ! 6: Sand Clay Loam
     ZOBLER_MEDIUM_FINE,    & ! 7: Clay Loam
     ZOBLER_MEDIUM,         & ! 8: Silty Clay Loam
     ZOBLER_COARSE_FINE,    & ! 9: Sandy Clay
     ZOBLER_FINE,           & !10: Silty Clay
     ZOBLER_FINE,           & !11: Clay
     ZOBLER_ORGANIC         & !12: Organic
  /)

contains

  !===========================================================================
  ! Subroutine: map_cptec_surface_types
  !===========================================================================
  !BOP
  ! !ROUTINE: map_cptec_surface_types
  ! !DESCRIPTION:
  !    Allocates and assigns the mapping arrays that convert model-specific
  !    vegetation classes into CRTM surface type codes for microwave (MWAVE)
  !    and infrared (IR) channels. The mapping depends on the chosen land
  !    surface model (SSiB or IBIS) and the CRTM IR classification scheme.
  ! !ARGUMENTS:
  !    idlsm - INTEGER, intent(in)
  !       Identifier of the land surface model:
  !         1 = SSiB, 2 = IBIS.
  ! !NOTES:
  !    - Uses global parameter arrays (ssib_to_gfs, ssib_to_npoess, etc.).
  !    - CRTM_IRlandCoeff_Classification() returns the IR scheme name
  !      ('NPOESS' or 'USGS') to select the appropriate IR mapping.
  ! !ERROR HANDLING:
  !    Prints an error and stops execution if idlsm is not 1 or 2.
  !EOP
  subroutine map_cptec_surface_types(idlsm)
    integer, intent(in) :: idlsm
    ! Mapping land surface type to CRTM surface fields

    if (idlsm == LSM_SSIB) then
       allocate(map_to_crtm_soil(SSIB_SOIL_N))
       map_to_crtm_soil = ssib_to_zobler

       allocate(map_to_crtm_mw(SSIB_VEG_N))
       map_to_crtm_mw = ssib_to_gfs

       allocate(map_to_crtm_ir(SSIB_VEG_N))
       select case (trim(CRTM_IRlandCoeff_Classification()))
          case ('NPOESS'); map_to_crtm_ir = ssib_to_npoess
          case ('IGBP');   map_to_crtm_ir = ssib_to_igbp
          case ('USGS');   map_to_crtm_ir = ssib_to_usgs
          case default; map_to_crtm_ir = ssib_to_gfs
       end select

    else if (idlsm == LSM_IBIS) then
       allocate(map_to_crtm_soil(IBIS_SOIL_N))
       map_to_crtm_soil = ibis_to_zobler

       allocate(map_to_crtm_mw(IBIS_VEG_N))
       map_to_crtm_mw = ibis_to_gfs

       allocate(map_to_crtm_ir(IBIS_VEG_N))
       select case (trim(CRTM_IRlandCoeff_Classification()))
          case ('NPOESS'); map_to_crtm_ir = ibis_to_npoess
          case ('IGBP');   map_to_crtm_ir = ibis_to_igbp
          case ('USGS');   map_to_crtm_ir = ibis_to_usgs
          case default; map_to_crtm_ir = ibis_to_gfs
       end select

    else
       write(6,*) "ERROR: Invalid land surface model. Only SSiB (1) or IBIS (2) are valid. idlsm=", idlsm
       call stop2(72)
    endif

  end subroutine map_cptec_surface_types


  !-------------------------------------------------------------------------------
  !BOP
  ! !ROUTINE: map_bam_to_crtm
  ! !DESCRIPTION:
  !    Maps vegetation, soil, water, and ice indices from the BAM model (SSiB/IBIS)
  !    into the CRTM_Surface_type structure, filling the corresponding fields:
  !    land, vegetation, soil, water, and ice.
  ! !INTERFACE:
  !    subroutine map_bam_to_crtm(itype, istype, sensor_type, sfc)
  ! !ARGUMENTS:
  !    itype       - integer, intent(in)
  !                  BAM index for vegetation/water/ice.
  !    istype      - integer, intent(in)
  !                  BAM index for soil type.
  !    sensor_type - integer, intent(in)
  !                  Sensor type used to identify microwave, IR, VIS, UV, etc.
  !    lai_type    - leaf-area-index for various types
  !    sfc         - type(CRTM_Surface_type), intent(inout)
  !                  CRTM structure to be populated with the mapped values.
  ! !REVISION HISTORY:
  !  08 Jul 2025 - J. G. de Mattos - Initial version.
  !EOP
  !BOC
  subroutine map_bam_to_crtm(itype, istype, sensor_type, lai_type, sfc)
    implicit none
    integer(i_kind), intent(in)            :: itype, istype
    integer(i_kind), intent(in)            :: sensor_type
    integer(i_kind), intent(inout)         :: lai_type
    type(CRTM_Surface_type), intent(inout) :: sfc
  
    integer(i_kind) :: ir_idx, mw_idx, soil_idx
    logical :: ok_ir, ok_mw, ok_soil
    integer, parameter :: VT = WATER_TYPE, IT = ICE_TYPE
  
    ! 1) Initialize all subfields to default
    sfc%Land_Type       = DEFAULT_LAND_TYPE
    sfc%Vegetation_Type = DEFAULT_VEGETATION_TYPE
    sfc%Soil_Type       = DEFAULT_SOIL_TYPE
    sfc%Water_Type      = DEFAULT_WATER_TYPE
    sfc%Ice_Type        = DEFAULT_ICE_TYPE
  
    ! 2) Compute CRTM indices
    ir_idx   = map_to_crtm_ir(safe_index(itype,map_to_crtm_ir))
    mw_idx   = map_to_crtm_mw(safe_index(itype,map_to_crtm_mw))
    soil_idx = map_to_crtm_soil(safe_index(istype, map_to_crtm_soil))

    ! 3) Validate index ranges
    ok_ir   = is_valid_veg(ir_idx)
    ok_mw   = is_valid_veg(mw_idx)
    ok_soil = is_valid_soil(soil_idx)
 
    ! 4) Special cases (water and ice) — early return
    if (ir_idx == VT) then
      sfc%Water_Type = 1  ! sea water (IR/VIS)
      return
    else if (ir_idx == IT) then
      sfc%Ice_Type = 1    ! new ice (IR/VIS)
      return
    end if
  
    ! 5) Invalid land case — fallback to compacted soil
    if (.not.(ok_ir .and. ok_mw .and. ok_soil)) then
      sfc%Land_Type       = DEFAULT_LAND_TYPE
      sfc%Soil_Type       = DEFAULT_SOIL_TYPE
      sfc%Vegetation_Type = DEFAULT_VEGETATION_TYPE
      return
    end if
  
    ! 6) Valid land case — fill land/veg/soil
    if (sensor_type == microwave_sensor) then
      sfc%Land_Type = ir_idx
    else
      sfc%Land_Type = mw_idx
    end if
    sfc%Vegetation_Type = mw_idx
    sfc%Soil_Type       = soil_idx
    lai_type            = mw_idx

    if (sfc%land_type < 0 ) then
      write(6,*) 'ERROR: land_type fora do intervalo [0,', max(SSIB_VEG_N,IBIS_VEG_N),'] =', sfc%land_type
      stop 'Invalid CRTM mapping'
    endif
  end subroutine map_bam_to_crtm
  !EOC
  
  !-------------------------------------------------------------------------------
  !BOP
  ! !ROUTINE: is_valid_veg
  ! !DESCRIPTION:
  !    Checks whether the vegetation index is within the valid range
  !    defined by either the SSiB or IBIS schemes.
  ! !INTERFACE:
  !    pure function is_valid_veg(idx) result(ok)
  ! !ARGUMENTS:
  !    idx - integer, intent(in)
  !          Vegetation index to be validated.
  ! !RETURN VALUE:
  !    ok  - logical
  !          .true. if idx is within [1, SSIB_VEG_N] or [1, IBIS_VEG_N];
  !          .false. otherwise.
  ! !REVISION HISTORY:
  !  08 Jul 2025 - J. G. de Mattos - Initial version.
  !EOP
  pure function is_valid_veg(idx) result(ok)
    implicit none
    integer, intent(in) :: idx
    logical :: ok
    ok = (idx >= 1 .and. idx <= SSIB_VEG_N) .or. &
         (idx >= 1 .and. idx <= IBIS_VEG_N)
  end function is_valid_veg
  
  !-------------------------------------------------------------------------------
  !BOP
  ! !ROUTINE: is_valid_soil
  ! !DESCRIPTION:
  !    Checks whether the soil index is within the valid range
  !    defined by either the SSiB or IBIS schemes.
  ! !INTERFACE:
  !    pure function is_valid_soil(idx) result(ok)
  ! !ARGUMENTS:
  !    idx - integer, intent(in)
  !          Soil index to be validated.
  ! !RETURN VALUE:
  !    ok  - logical
  !          .true. if idx is within [1, SSIB_SOIL_N] or [1, IBIS_SOIL_N];
  !          .false. otherwise.
  ! !REVISION HISTORY:
  !  08 Jul 2025 - J. G. de Mattos - Initial version.
  !EOP
  pure function is_valid_soil(idx) result(ok)
    implicit none
    integer, intent(in) :: idx
    logical :: ok
    ok = (idx >= 1 .and. idx <= SSIB_SOIL_N) .or. &
         (idx >= 1 .and. idx <= IBIS_SOIL_N)
  end function is_valid_soil

  !BOP
  ! !ROUTINE: safe_index
  ! !DESCRIPTION:
  !   This function ensures safe indexing into a dynamically allocated Fortran array.
  !   It returns a valid index (`clipped_idx`) that lies within the bounds of the input array.
  !   If the provided index `idx` is smaller than the lower bound of the array,
  !   the lower bound is returned. If it is greater than the upper bound, the upper bound is returned.
  !   Otherwise, the index is returned unchanged.
  !
  !   This is particularly useful when working with arrays that are allocated at runtime,
  !   and whose bounds may not be known in advance.
  !
  ! !INTERFACE:
  !   integer function safe_index(array, idx) result(clipped_idx)
  !
  ! !ARGUMENTS:
  !   integer,           intent(in)           :: idx
  !     → Index value that may be out of bounds.
  !   integer, dimension(:), intent(in)       :: array
  !     → One-dimensional array with dynamic bounds.
  !
  ! !RETURN VALUE:
  !   integer :: clipped_idx
  !     → Index guaranteed to be within the range [lbound(array,1), ubound(array,1)].
  !
  ! !IMPLEMENTATION DETAILS:
  !   - Uses the intrinsic Fortran functions `lbound` and `ubound` to retrieve
  !     the lower and upper bounds of the array.
  !   - Combines `max` and `min` intrinsics to clip the index efficiently:
  !       `clipped_idx = max(imin, min(idx, imax))`
  !
  ! !EXAMPLE USAGE:
  !   integer, allocatable :: a(:)
  !   integer :: i
  !   allocate(a(3:10))
  !   i = safe_index(a, 1)   ! returns 3
  !   i = safe_index(a, 11)  ! returns 10
  !   i = safe_index(a, 5)   ! returns 5
  !
  ! !REVISION HISTORY:
  !   09 Jul 2025 - J. G. de Mattos - Initial version.
  !EOP
  function safe_index(idx, array) result(clipped_idx)
    implicit none
    integer, intent(in)           :: idx
    integer, dimension(:), intent(in) :: array
    integer :: clipped_idx
    integer :: imin, imax
  
    imin = lbound(array, 1)
    imax = ubound(array, 1)
  
    clipped_idx = max(imin, min(idx, imax))
  end function safe_index
  !-------------------------------------------------------------------------------
  !BOP
  ! !ROUTINE: finalize_cptec_crtm_mapping
  ! !DESCRIPTION:
  !   This subroutine deallocates the global mapping arrays used for translating
  !   CPTEC/BAM vegetation and soil types into CRTM surface classification schemes.
  !   It should be called once the mappings are no longer needed (e.g., at the end
  !   of a model initialization or processing step) to free memory.
  !
  ! !USAGE:
  !   call finalize_cptec_crtm_mapping()
  !
  ! !SIDE EFFECTS:
  !   The global arrays `map_to_crtm_ir`, `map_to_crtm_mw`, and `map_to_crtm_soil`
  !   will be deallocated if currently allocated.
  !
  ! !REVISION HISTORY:
  !   09 Jul 2025 - J. G. de Mattos - Initial version.
  !EOP
  subroutine finalize_cptec_crtm_mapping()
    implicit none
  
    if (allocated(map_to_crtm_ir)) then
      deallocate(map_to_crtm_ir)
    end if
  
    if (allocated(map_to_crtm_mw)) then
      deallocate(map_to_crtm_mw)
    end if
  
    if (allocated(map_to_crtm_soil)) then
      deallocate(map_to_crtm_soil)
    end if
  end subroutine finalize_cptec_crtm_mapping

end module cptecbam_crtm_mapping

