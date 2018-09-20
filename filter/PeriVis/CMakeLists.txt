
option(BUILD_BONDSWRITER "BondsWriter bauen" ON)
IF( BUILD_BONDSWRITER )
  add_subdirectory( BondsWriter)
ENDIF()

option(BUILD_TETRIDGEFILTER "TetRidgeFilter bauen" ON)
IF( BUILD_TETRIDGEFILTER )
  add_subdirectory( TetRidgeFilter)
ENDIF()

option(BUILD_GENERATEBONDS "GenerateBonds bauen" ON)
IF( BUILD_GENERATEBONDS )
  add_subdirectory( GenerateBonds)
ENDIF()

option(BUILD_POLYPOINTSET "PolyPointSet bauen" ON)
IF( BUILD_POLYPOINTSET )
  add_subdirectory( PolyPointSet)
ENDIF()

option(BUILD_CLIPLINES "ClipLines bauen" ON)
IF( BUILD_CLIPLINES )
  add_subdirectory( ClipLines)
ENDIF()

option(BUILD_LAMMPSTEMPORALREADER "LAMMPSTemporalReader bauen" ON)
IF( BUILD_LAMMPSTEMPORALREADER )
  add_subdirectory( LAMMPSTemporalReader)
ENDIF()

option(BUILD_LAMMPSTEMPORALREADER_MF "LAMMPSTemporalReaderMultiFit bauen" OFF)
IF( BUILD_LAMMPSTEMPORALREADER_MF )
  add_subdirectory( LAMMPSTemporalReaderMultiFit )
ENDIF()

option(BUILD_LAMMPSWRITER "LAMMPSWriter bauen" ON)
IF( BUILD_LAMMPSWRITER )
  add_subdirectory( LAMMPSWriter )
ENDIF()

option(BUILD_LIGHTSCATTERING "LightScattering bauen" ON)
IF( BUILD_LIGHTSCATTERING )
  add_subdirectory( LightScattering )
ENDIF()

option(BUILD_PERIPCA "PeriPCA bauen" ON)
IF( BUILD_PERIPCA )
  add_subdirectory( PeriPCA )
ENDIF()

option(BUILD_PERIMULTIFIT "PeriMultiFit bauen" OFF)
IF( BUILD_PERIMULTIFIT )
  add_subdirectory( PeriMultiFit )
ENDIF()

option(BUILD_PERIMULTIFIT_TAYLOR "PeriMultiFitTaylor bauen" OFF)
IF( BUILD_PERIMULTIFIT_TAYLOR )
  add_subdirectory( PeriMultiFitTaylor )
ENDIF()

option(BUILD_PERIMULTIFIT_HESSE "PeriMultiFitHesse bauen" OFF)
IF( BUILD_PERIMULTIFIT_HESSE )
  add_subdirectory( PeriMultiFitHesse )
ENDIF()

option( BUILD_CLIPPOLYDATA "ClipPolyData bauen" ON)
IF( BUILD_CLIPPOLYDATA )
  add_subdirectory( ClipPolyData )
ENDIF()

option(BUILD_IMPLICITPOLYDATACLIP "ImplicitPolydataClip bauen" ON)
IF( BUILD_IMPLICITPOLYDATACLIP )
  add_subdirectory( ImplicitPolydataClip )
ENDIF()

option(BUILD_NEARESTNEIGHBORSAMPLING "NearestNeighborSampling bauen" OFF)
IF( BUILD_NEARESTNEIGHBORSAMPLING )
  add_subdirectory( NearestNeighborSampling)
ENDIF()

option(BUILD_COUNTBONDS "CountBonds bauen" ON)
IF( BUILD_COUNTBONDS )
  add_subdirectory( CountBonds)
ENDIF()

option(BUILD_POWERCRUSTSURFACERECONSTRUCTION "PowerCrustSurfaceReconstruction bauen" OFF)
IF( BUILD_POWERCRUSTSURFACERECONSTRUCTION )
  add_subdirectory( PowerCrustSurfaceReconstruction )
ENDIF()

option(BUILD_CONVERTTOPOLYDATA "ConvertToPolydata bauen" ON)
IF( BUILD_CONVERTTOPOLYDATA )
  add_subdirectory( ConvertToPolydata )
ENDIF()

option(BUILD_CLEANBONDS "CleanBonds bauen" ON)
IF( BUILD_CLEANBONDS )
  add_subdirectory( CleanBonds )
ENDIF()

option(BUILD_PERICONNECTEDCOMPONENTS "PeriConnectedComponents bauen" ON)
IF( BUILD_PERICONNECTEDCOMPONENTS )
  add_subdirectory( PeriConnectedComponents )
ENDIF()

option(BUILD_PERISTRESS "PeriStress bauen" ON)
IF( BUILD_PERISTRESS )
  add_subdirectory( PeriStress )
ENDIF()

option(BUILD_MESHDIFFERENCE "MeshDifference bauen" ON)
IF( BUILD_MESHDIFFERENCE )
  add_subdirectory( MeshDifference )
ENDIF()

option(BUILD_LOADBONDS "LoadBonds bauen" ON)
IF( BUILD_LOADBONDS )
  add_subdirectory( LoadBonds )
ENDIF()

option(BUILD_EXTRACTBROKENBONDS "ExtractBrokenBonds bauen" ON)
IF( BUILD_EXTRACTBROKENBONDS )
  add_subdirectory( ExtractBrokenBonds )
ENDIF()

option(BUILD_TETRIDGE "TetRidge bauen" ON)
IF( BUILD_TETRIDGE )
  add_subdirectory( TetRidge )
ENDIF()

option(BUILD_TETRIDGE_2D "TetRidge2D bauen" ON)
IF( BUILD_TETRIDGE_2D )
  add_subdirectory( TetRidge2D )
ENDIF()

option(BUILD_TETSMOOTH "TetSmooth bauen" ON)
IF( BUILD_TETSMOOTH )
  add_subdirectory( TetSmooth )
ENDIF()

option(BUILD_TETGEN "TetGen bauen" ON)
IF( BUILD_TETGEN )
  add_subdirectory( TetGen )
ENDIF()

option(BUILD_TETTRIMSURFACE "TetTrimSurface bauen" ON)
IF( BUILD_TETTRIMSURFACE )
  add_subdirectory( TetTrimSurface )
ENDIF()

option(BUILD_TETCLIP "TetClip bauen" ON)
IF( BUILD_TETCLIP )
  add_subdirectory( TetClip )
ENDIF()

option(BUILD_IMPLICITTETCLIP "ImplicitTetClip bauen" ON)
IF( BUILD_IMPLICITTETCLIP )
  add_subdirectory( ImplicitTetClip )
ENDIF()

#add_subdirectory(ConvertToPolydata)