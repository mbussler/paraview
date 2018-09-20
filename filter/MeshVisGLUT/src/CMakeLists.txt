set( header
	BoundingBox.h
	Camera.h
	DisplayFPS.h
	InlineFunctions.h
	KdTree.h
	KdTreeCuda.h
	KdTreeTest.h
	Matrix.h
	MeshLoader.h
	MeshProcessor.h
	MeshVis.h
	Particle.h
	ParticleImplementationCuda.h
	ParticleSeeder.h
	ParticlesGL.h
	quickselect.h
	RecordingTools.h
	shaders.h
	TetMesh.h
	TetMeshCuda.h
	TetMeshGL.h
	Tetrahedron.h
	TetrahedronCreator.h
	TetWalker.h
	timer.h
	TimeStepControllerCuda.h
	TimestepCuda.h
	Typedefs.h
	Vector.h
	VTKLoader.h
)

set( sources
	BoundingBox.cpp
	Camera.cpp
	KdTree.cpp
	KdTreeCuda.cpp
	KdTreeTest.cpp
	MeshLoader.cpp
	MeshProcessor.cpp
	MeshVis.cpp
	Particle.cpp
	ParticleImplementationCuda.cpp
	ParticlesGL.cpp
	RecordingTools.cpp
	shaders.cpp
	TetMesh.cpp
	TetMeshCuda.cpp
	TetMeshGL.cpp
	Tetrahedron.cpp
	TetrahedronCreator.cpp
	TetWalker.cpp
	TimeStepControllerCuda.cpp
	TimestepCuda.cpp
	visMain.cpp
	VTKLoader.cpp
)

set( header_cuda 
	common.cuh
	KdTraverse.cuh
	Particles.cuh
	Tetwalk.cuh
)

set( source_cuda
	common.cu
	KdTraverse.cu
	MeshFunctions.cu
	Particles.cu
	Tetwalk.cu
)

set( kernels_cuda
	KdTraverse_kernel.cu
	Integration_kernel.cu
	Particles_kernel.cu
	TetWalk_kernel.cu
)

# main.cpp
# mainCPU.cpp
# main_test.cpp

#project(MeshVisGLUT)

IF( USE_TETGEN )
	set( sources ${sources} TetgenLoader.cpp )
	set( header  ${header}  TetgenLoader.h )
ENDIF()

add_definitions( -DNOMINMAX )
CUDA_ADD_EXECUTABLE( MeshVisGLUT ${sources} ${source_cuda} )
target_link_libraries( MeshVisGLUT ${OPENGL_LIBRARIES} ${VTK_LIBRARIES} ${Boost_LIBRARIES} ${GLEW_LIBRARIES} ${GLUT_LIBRARIES} ${CUDA_CUDA_LIBRARY})

INSTALL( TARGETS MeshVisGLUT  RUNTIME DESTINATION bin
                            LIBRARY DESTINATION bin)

IF( USE_TETGEN )
	target_link_libraries( MeshVisGLUT optimized tetgen debug tetgen_d )
ENDIF()
