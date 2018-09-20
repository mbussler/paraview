#include "TimestepCPU.h"

TimestepCPU::TimestepCPU()
{
	pointInTime = 0.0f;
	timestepIndex=0;

	m_meshgl.reset( new TetMeshGL() );
};

TimestepCPU::~TimestepCPU()
{
};

void TimestepCPU::init(int _timestepIndex)
{
	timestepIndex = _timestepIndex;

	printf("\n");
	printf("Initializing [Timestep %d]:\n", timestepIndex);
	printf("Initializing TetrahedronCreator..");
		m_tetCreator.reset( new TetrahedronCreator() );
		m_tetCreator->setMesh( m_meshgl );
		m_tetCreator->init();
	printf("done.\n");

	// lazy evaluation..
	//printf("Initializing Vertex-Buffer Objects..", timestepIndex);
	//	m_meshgl->initVBOs();
	//printf("done.\n");
	m_meshgl->init();

	printf("Processing Kd-Tree..\n", timestepIndex);
		m_kdtree.reset( new KdTreeTest() );
		m_kdtree->createKdTree( m_meshgl );

	m_tetWalker.reset( new TetWalker() );
	m_tetWalker->setTetrahedronCreator( m_tetCreator );
	m_tetWalker->setKDTree( m_kdtree );

};

void TimestepCPU::locateParticles( ParticleList& pl )
{
	ParticleListIterator pli;
	int cellIndex;
	for( pli=pl.begin(); pli!=pl.end(); pli++)
	{
		// dismiss particles that are out of the field
		if( pli->bOutOfField )
			continue;

		// use the tetwalker to locate the cell for the given position
		cellIndex = m_tetWalker->locateCell( pli->Position );

		// store the cell index in the particles Cells array
		pli->Cell[ timestepIndex ] = cellIndex;

		//m_meshgl.getStartCellsList().push_back(cellIndex);
	};
};

void TimestepCPU::locateParticle( Particle& p )
{
	// use the tetwalker to locate the cell for the given position
	p.Cell[ timestepIndex ] = m_tetWalker->locateCell( p.Position );
};

void TimestepCPU::interpolateVelocityAt( const Point& position, int& startCell, Vector& velocity)
{
	// store starting index in MeshGL's starting cells list
	// m_meshgl.getStartCellsList().push_back(startCell);

	// calculate interpolated Velocity for the position starting the TetWalk at the startCell
	m_tetWalker->interpolateVelocityAt( position, startCell, velocity );
};

void TimestepCPU::storeTraversedCells( bool store )
{
	if( store )
	{
		m_tetWalker->setupVis( m_meshgl->getTraversedCellsList() );
	}
	else
	{
		m_meshgl->resetTraversedCells();
		m_tetWalker->setupVis( 0 );
	}
}
