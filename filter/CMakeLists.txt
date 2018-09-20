cmake_minimum_required(VERSION 2.8)

PROJECT(PV-Plugins)

option(BUILD_GROWTHVELOCITY "GrowthVelocity bauen" ON)
IF( BUILD_GROWTHVELOCITY )
  add_subdirectory( GrowthVelocity)
ENDIF()

option(BUILD_RIDGEANGLE "RidgeAngle bauen" ON)
IF( BUILD_RIDGEANGLE )
  add_subdirectory( RidgeAngle)
ENDIF()

option(BUILD_CELLSIZE "CellSize bauen" ON)
IF( BUILD_CELLSIZE )
  add_subdirectory( CellSize)
ENDIF()

option(BUILD_WRITEPOLYDATA "WritePolyData bauen" ON)
IF( BUILD_WRITEPOLYDATA )
  add_subdirectory( WritePolyData)
ENDIF()

option(BUILD_TIMEGEOMETRYUNSTRUCTURED "TimeGeometryUnstructured bauen" ON)
IF( BUILD_TIMEGEOMETRYUNSTRUCTURED )
  add_subdirectory( TimeGeometryUnstructured)
ENDIF()

option(BUILD_TIMEGEOMETRY "TimeGeometry bauen" ON)
IF( BUILD_TIMEGEOMETRY )
  add_subdirectory( TimeGeometry)
ENDIF()

option(BUILD_ADDTIME "AddTime bauen" ON)
IF( BUILD_ADDTIME )
  add_subdirectory( AddTime)
ENDIF()

option(BUILD_VTKRIDGE "vtkRidge bauen" ON)
IF( BUILD_VTKRIDGE )
  add_subdirectory( vtkRidge)
ENDIF()

option(BUILD_MESHCOMPONENTS "MeshComponents bauen" ON)
IF( BUILD_MESHCOMPONENTS )
  add_subdirectory( MeshComponents)
ENDIF()

option(BUILD_GAUSSIANSMOOTH "GaussianSmooth bauen" ON)
IF( BUILD_GAUSSIANSMOOTH )
  add_subdirectory( GaussianSmooth)
ENDIF()

option(BUILD_RIDGECOMPARE "RidgeCompare bauen" ON)
IF( BUILD_RIDGECOMPARE )
  add_subdirectory( RidgeCompare )
ENDIF()

option(BUILD_RIDGEMERGE "RidgeMerge bauen" ON)
IF( BUILD_RIDGEMERGE )
  add_subdirectory( RidgeMerge )
ENDIF()

option(BUILD_RIDGEGROW "RidgeGrow bauen" ON)
IF( BUILD_RIDGEGROW )
  add_subdirectory( RidgeGrow )
ENDIF()

option(BUILD_RIDGEGROWSTRESS "RidgeGrowStress bauen" ON)
IF( BUILD_RIDGEGROWSTRESS )
  add_subdirectory( RidgeGrowStress )
ENDIF()

option(BUILD_RIDGESTEADY "RidgeSteady bauen" ON)
IF( BUILD_RIDGESTEADY )
  add_subdirectory( RidgeSteady )
ENDIF()

option(BUILD_RIDGETIME "RidgeTime bauen" ON)
IF( BUILD_RIDGETIME )
  add_subdirectory( RidgeTime )
ENDIF()

option(BUILD_LEASTSQUARESGRADIENTS "LeastSquaresGradients bauen" ON)
IF( BUILD_LEASTSQUARESGRADIENTS )
  add_subdirectory( LeastSquaresGradients )
ENDIF()

option(BUILD_DATRAWTEMPORALREADER "DatRawTemporalReader bauen" OFF)
IF( BUILD_DATRAWTEMPORALREADER )
  add_subdirectory( DatRawTemporalReader)
ENDIF()

option(BUILD_PERIVIS_PLUGINS "PeriVis Plugins bauen" ON)
IF( BUILD_PERIVIS_PLUGINS )
  add_subdirectory( PeriVis )
ENDIF()

option(BUILD_CELLFLUX "cellFlux bauen" ON)
IF(BUILD_CELLFLUX ) 
  add_subdirectory( cellFlux )
ENDIF()

option(BUILD_CUSTOMDATASET "customDataSet bauen" OFF)
IF(BUILD_CUSTOMDATASET ) 
  add_subdirectory( customDataSet )
ENDIF()

option(BUILD_DIVERGENCE "divergence bauen" ON)
IF(BUILD_DIVERGENCE ) 
  add_subdirectory( divergence )
ENDIF()

option(BUILD_EIGENANALYSIS "eigenanalysis bauen" OFF)
IF(BUILD_EIGENANALYSIS ) 
  add_subdirectory( eigenanalysis )
ENDIF()

option(BUILD_EIGENVALUETEST "EigenValueTest bauen" OFF)
IF(BUILD_EIGENVALUETEST ) 
  add_subdirectory( EigenValueTest )
ENDIF()

option(BUILD_EIGENVECTORS "eigenvectors bauen" ON)
IF(BUILD_EIGENVECTORS ) 
  add_subdirectory( eigenvectors )
ENDIF()

option(BUILD_MAGNITUDE "magnitude bauen" ON)
IF(BUILD_MAGNITUDE ) 
  add_subdirectory( magnitude )
ENDIF()

option(BUILD_MICHELLEVY "MichelLevy bauen" ON)
IF(BUILD_MICHELLEVY ) 
  add_subdirectory( MichelLevy )
ENDIF()

option(BUILD_RATEOFSTRAIN "rateOfStrain bauen" ON)
IF(BUILD_RATEOFSTRAIN )
  add_subdirectory( rateOfStrain )
ENDIF()

option(BUILD_SHEPARDMETHOD "ShepardMethod bauen" ON)
IF(BUILD_SHEPARDMETHOD )
  add_subdirectory( ShepardMethod )
ENDIF()

option(BUILD_STREAMTRACEREV "StreamTracerEV bauen" OFF)
IF(BUILD_STREAMTRACEREV )
  add_subdirectory( StreamTracerEV )
ENDIF()

option(BUILD_TENSORNORM "TensorNorm bauen" ON)
IF(BUILD_TENSORNORM )
  add_subdirectory( TensorNorm )
ENDIF()

option(BUILD_VECTOR2DTO3D "vector2Dto3D bauen" ON)
IF(BUILD_VECTOR2DTO3D )
  add_subdirectory( vector2Dto3D )
ENDIF()

option(BUILD_VECTORPCA "vectorPCA bauen" ON)
IF(BUILD_VECTORPCA )
  add_subdirectory( vectorPCA )
ENDIF()

option(BUILD_STRAINRATETENSORGLYPHS "StrainRateTensorGlyphs bauen" ON)
IF(BUILD_STRAINRATETENSORGLYPHS )
  add_subdirectory( StrainRateTensorGlyphs )
ENDIF()

option(BUILD_STREAMTRACERPCA "StreamTracerPCA bauen" OFF)
IF(BUILD_STREAMTRACERPCA ) 
  add_subdirectory( StreamTracerPCA )
ENDIF()

option(BUILD_TENSORWRITERAW "TensorWriteRaw bauen" ON)
IF(BUILD_TENSORWRITERAW )
  add_subdirectory( TensorWriteRaw )
ENDIF()

option(BUILD_VECTORCALC "vectorCalc bauen" ON)
IF(BUILD_VECTORCALC )
  add_subdirectory( vectorCalc )
ENDIF()

option(BUILD_WARP2D "warp2D bauen" ON)
IF(BUILD_WARP2D )
  add_subdirectory( warp2D )
ENDIF()

option(BUILD_MAKETENSOR "MakeTensor bauen" ON)
IF(BUILD_MAKETENSOR )
  add_subdirectory( MakeTensor )
ENDIF()

option(BUILD_VONMIESESSTRESS "VonMiesesStress bauen" ON)
IF(BUILD_VONMIESESSTRESS )
  add_subdirectory( VonMiesesStress )
ENDIF()

option(BUILD_POWERCRUST "PowerCrust bauen" OFF)
IF(BUILD_POWERCRUST )
  add_subdirectory( PowerCrust )
ENDIF()

option(BUILD_TENSORGLYPHSNORMALIZED "TensorGlyphsNormalized bauen" ON)
IF(BUILD_TENSORGLYPHSNORMALIZED )
  add_subdirectory( TensorGlyphsNormalized )
ENDIF()
