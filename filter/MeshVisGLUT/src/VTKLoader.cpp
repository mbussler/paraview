#include "VTKLoader.h"
#include "InlineFunctions.h"

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkDataSetAttributes.h>
#include <vtkCleanPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkQuadraticTetra.h>
#include <vtkTetra.h>

/// constructor
VTKLoader::VTKLoader() :
	a(0.0f), b(1.0f), c(0.1f),
	m_domain( BoundingBox( Point( -1, -1, -1), Point( 1, 1, 1)))
{
	
	
};

/// destructor
VTKLoader::~VTKLoader()
{
};

/// load nodes from file
void VTKLoader::loadFromFile( const char* filename, TetMeshPtr mesh)
{
	sprintf( mesh->filename, "%s", filename);

	// check if file exists
	ifstream file(filename);
	if(!file)
	{
		printf("[VTKLoader] Error opening %s: File not found!\n", filename);
		exit(0);
	}
	else
	{
		file.close();
	}

	vtkSmartPointer<vtkUnstructuredGrid> grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
		vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();

	printf("\n");
	printf("[VTKLoader] Loading Mesh from %s.\n", filename);

	reader->SetFileName( filename );
	reader->Update();
	grid = reader->GetOutput();

    grid->GetPointData()->SetActiveVectors("wf");

	VtkGridToMesh( grid, mesh);

	if( true || !loadNeighborsFromFile(filename, mesh) )
	{
		calculateNeighbors( grid, mesh );
		saveNeighborsToFile( filename, mesh );
	}

};

void VTKLoader::calculateNeighbors( vtkUnstructuredGrid *grid, TetMeshPtr mesh)
{
	assert( mesh->cells_count > 0);

	printf("[VTKLoader] Calculating neighbors");

	mesh->neighbors = new Cell[ mesh->cells_count ];

    const int print_each = (mesh->cells_count / 10)+1;

	list<Face> outer_faces;
	list<Vector> normals;
    list<Vector> face_vertices;

    vtkIdType cellId, faceId;
    
    for( cellId = 0; cellId < mesh->cells_count; cellId++)
    {
      vtkCell* cell = grid->GetCell(cellId);
      
      for( faceId = 0; faceId<cell->GetNumberOfFaces(); faceId++)
      {
        
        vtkCell* face = cell->GetFace(faceId);

        vtkSmartPointer<vtkIdList> facePointIds = vtkSmartPointer<vtkIdList>::New();
        for( vtkIdType ptId = 0; ptId < face->GetNumberOfPoints(); ptId++) {
          facePointIds->InsertNextId( face->GetPointId(ptId));
        }

        vtkSmartPointer<vtkIdList> neighborIds = vtkSmartPointer<vtkIdList>::New();
        grid->GetCellNeighbors(cellId, facePointIds, neighborIds );
        
        int neighbor = -1;

        int numNeighborIds = neighborIds->GetNumberOfIds() ;
        if( neighborIds->GetNumberOfIds() > 0)
        {
            neighbor = neighborIds->GetId(0);
        }
        else
        {
            Face f;
            f.pts[0] = facePointIds->GetId(0);
            f.pts[1] = facePointIds->GetId(1);
            f.pts[2] = facePointIds->GetId(2);              

            // calculate face normal
            Vector p1 = Pt2Vec( mesh->nodes[f.pts[0]] );
            Vector p2 = Pt2Vec( mesh->nodes[f.pts[1]] );
            Vector p3 = Pt2Vec( mesh->nodes[f.pts[2]] );

            Vector f_center = 1.0f/3.0f * ( p1 + p2 + p3);
            Vector normal = VectorNormalize(VectorCross( p1-p2, p2-p3));

            // calculate cell center
            Cell c = mesh->cells[cellId];
            Vector vs[4];

            for( int j=0; j<4; j++) {
                vs[j] = Pt2Vec( mesh->nodes[ mesh->cells[cellId].indices[j] ] );
            }

            Vector c_center = 1.0f/4.0f * (vs[0]+vs[1]+vs[2]+vs[3]);

            if( VectorDot( (c_center - f_center), normal) < 0 )
            {
                normal = normal * -1;
                int swap = f.pts[1];
                f.pts[1] = f.pts[2];
                f.pts[2] = swap;
                
                Vector swp = p1;
                p1 = p2; 
                p2 = swp;
            }

            face_vertices.push_back(p1);
            face_vertices.push_back(p2);
            face_vertices.push_back(p3);
            
            outer_faces.push_back(f);
            normals.push_back(normal);
        }
        mesh->neighbors[cellId].indices[faceId] = neighbor;
      }

      if( cellId%print_each == 0 ) printf(".");
	}

	mesh->copyOuterFaces(outer_faces);
	mesh->copyNormals(normals);
    mesh->copyFaceVertices(face_vertices);
    
    printf("done.\n");

};

void VTKLoader::createRandomPointSet( TetMeshPtr mesh, int num_points, bool save /*= true*/ )
{
	
	delete[] mesh->nodes;
	mesh->nodes = new Vertex[num_points];
	mesh->node_count = num_points;

	delete[] mesh->node_attributes;
	mesh->node_attributes = new VertexAttribute[num_points];
	mesh->attribute_count = num_points;

	mesh->updateBoundingBox( false );
	mesh->box = m_domain;
	mesh->grid = Point( 1, 1, 1 );

	printf("[VTKLoader] Create Random Pointset.\n");

	// Create the geometry of a point (the coordinate)
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();
	
	// create random point set and calculate flow field
	for( int i=0; i<num_points; i++)
	{
		// create random point in domain
		Point rnd = m_domain.RandomPointInBox();

		// store point in mesh
		mesh->nodes[i] = rnd.toFloat4();

		// store point in vtkPoints
		points->InsertNextPoint( rnd.c );

		// calculate flow at position
		Point fl = flow( rnd.x, rnd.y, rnd.z );
		mesh->node_attributes[i] = fl.toFloat4();
	}

	// Create a pointset object and add the points to it.
	vtkSmartPointer<vtkPointSet> pointset = 
		vtkSmartPointer<vtkPolyData>::New();
	pointset->SetPoints(points);

	printf("[VTKLoader] Calculate cells.\n");

	// Generate a tetrahedral mesh from the input points. By
	// default, the generated volume is the convex hull of the points.
	vtkSmartPointer<vtkDelaunay3D> delaunay3D =
		vtkSmartPointer<vtkDelaunay3D>::New();
	delaunay3D->SetInputData( pointset );
	delaunay3D->Update();

	char filename[50];
	sprintf(filename, "RandomPointSet_%d_a%4.2f_b%4.2f_c%4.2f.xml",num_points, a,b,c);

	if ( save )
	{
		// save mesh to file
		// Write the mesh as an unstructured grid
		vtkSmartPointer<vtkXMLDataSetWriter> writer =
			vtkSmartPointer<vtkXMLDataSetWriter>::New();
		writer->SetFileName ( filename );
		writer->SetInputConnection ( delaunay3D->GetOutputPort() );
		writer->Write();
	}

	vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = delaunay3D->GetOutput();

	// copy cells to mesh
	delete mesh->cells;
	int cell_count = vtkGrid->GetNumberOfCells();
	mesh->cells = new Cell[ cell_count ];
	mesh->cells_count = cell_count;

	// copy point ids
	for (int i=0; i<cell_count; i++ )
	{
		vtkIdList* cellids = vtkGrid->GetCell( i )->GetPointIds();
		for (int j = 0; j < 4; j++)
			mesh->cells[i].indices[j] = cellids->GetId(j);
	}

	calculateNeighbors( vtkGrid, mesh );

	if( save )
	{
		saveNeighborsToFile( filename, mesh );
	}

	printf("[VTKLoader] %d random nodes created.\n",num_points);

}

void VTKLoader::createTetrahedralGrid( TetMeshPtr mesh, int gridX, int gridY, int gridZ, float variance, bool save_to_file/*=false*/ )
{
	int num_points = gridX * gridY * gridZ;

	printf("[VTKLoader] Create Tetrahedral Grid.\n");

	// Create the geometry of a point (the coordinate)
	vtkSmartPointer<vtkPoints> points_in =
		vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkFloatArray> point_flow =
		vtkSmartPointer<vtkFloatArray>::New();

	point_flow->SetName("velocity");
	point_flow->SetNumberOfComponents(3);

	for( int k = 0; k < gridZ; k++)
	{
		for( int j = 0; j < gridY; j++)
		{
			for( int i = 0; i < gridX; i++)
			{
				//Point p = Point(i,j,k);
				Point p;
				p.x = ((REAL)i / (REAL)(gridX-1)) * m_domain.getLengthOfDim(0) + m_domain.point_min.x;
				p.y = ((REAL)j / (REAL)(gridY-1)) * m_domain.getLengthOfDim(1) + m_domain.point_min.y;
				p.z = ((REAL)k / (REAL)(gridZ-1)) * m_domain.getLengthOfDim(2) + m_domain.point_min.z;

				// store point in vtkPoints
				points_in->InsertNextPoint( p.c );

				// create random offset
				//p.x += ((i==0) || (i==gridX-1)) ? 0 : Random(-variance, variance);
				//p.y += ((j==0) || (j==gridY-1)) ? 0 : Random(-variance, variance);
				//p.z += ((k==0) || (k==gridZ-1)) ? 0 : Random(-variance, variance);

				// calculate synthetic flow
				Point f = flow( p.x, p.y, p.z);

				point_flow->InsertNextTuple3( f.x, f.y, f.z );
			}
		}
	}

	// Create a unstructured grid object and add the points to it.
	vtkSmartPointer<vtkUnstructuredGrid> vtkGrid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkGrid->SetPoints( points_in );
	vtkGrid->GetPointData()->AddArray(point_flow);

	printf("[VTKLoader] Calculate cells.\n");

	// Generate a tetrahedral mesh from the input points. By
	// default, the generated volume is the convex hull of the points.
	vtkSmartPointer<vtkDelaunay3D> delaunay3D =
		vtkSmartPointer<vtkDelaunay3D>::New();

	delaunay3D->SetInputData( vtkGrid );
	delaunay3D->Update();
	vtkGrid = delaunay3D->GetOutput();

	vtkCellArray *cells   = vtkGrid->GetCells();
	vtkPoints    *points  = vtkGrid->GetPoints();
	vtkPointData *data    = vtkGrid->GetPointData();
	vtkDataArray *da	  =	data->GetVectors("velocity");

	sprintf( mesh->filename, "TetMesh_%dx%dx%d_a%4.2f_b%4.2f_c%4.2f.vtu",gridX, gridY, gridZ, a,b,c);

	if ( save_to_file )
	{
		// save mesh to file
		// Write the mesh as an unstructured grid
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
			vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		writer->SetFileName ( mesh->filename );
		writer->SetInputData( vtkGrid );
		writer->Write();
	}

	mesh->updateBoundingBox( false );
	mesh->box = m_domain;
	mesh->grid = Point( gridX, gridY, gridZ );

	// copy to mesh
	VtkGridToMesh( vtkGrid, mesh);

	if( !loadNeighborsFromFile( mesh->filename, mesh) )
	{
		calculateNeighbors( vtkGrid, mesh );
		saveNeighborsToFile( mesh->filename, mesh );
	}

	printf("[VTKLoader] %d random nodes created.\n",num_points);
}

REAL VTKLoader::flow_u(REAL x, REAL y, REAL z)
{
	return a*x-b*y;
};

REAL VTKLoader::flow_v(REAL x, REAL y, REAL z)
{
	return b*x+a*y;
};

REAL VTKLoader::flow_w(REAL x, REAL y, REAL z)
{
	return -2*a*z+c;
};

Point VTKLoader::flow(REAL x, REAL y, REAL z)
{
	Point p;
	p.x = flow_u(x,y,z);
	p.y = flow_v(x,y,z);
	p.z = flow_w(x,y,z);
	return p;
};

void VTKLoader::VtkGridToMesh( vtkUnstructuredGrid* vtkGrid, TetMeshPtr mesh )
{
	vtkSmartPointer<vtkCellArray> cells = vtkGrid->GetCells();
	vtkSmartPointer<vtkPoints>    points  = vtkGrid->GetPoints();
	vtkSmartPointer<vtkPointData> data    = vtkGrid->GetPointData();

	mesh->node_count = points->GetNumberOfPoints();
	mesh->cells_count = cells->GetNumberOfCells();
	mesh->attribute_count = 3;

	printf("[VTKLoader] Loading %d nodes and %d cells with %d attribute(s) per node.\n", mesh->node_count, mesh->cells_count, mesh->attribute_count);
	printf("[VTKLoader] Loading mesh points..");

	delete[] mesh->nodes;
    int nNodes = (mesh->node_count / 10) * 4;
	mesh->nodes = new Vertex[ mesh->node_count ];

    float scale = 100.0f;
    
	for( int i=0; i< mesh->node_count; i++)
	{
		double* v = points->GetPoint(i);
		mesh->nodes[i].x = (float) v[0] * scale;
		mesh->nodes[i].y = (float) v[1] * scale;
		mesh->nodes[i].z = (float) v[2] * scale;
		mesh->nodes[i].w = 1.0f;

		mesh->box.updateBoundingBox(mesh->nodes[i]);
	}
	printf("done.\n");

	printf("[VTKLoader] Loading cells indices..");

	delete[] mesh->cells;
	mesh->cells     = new Cell[ mesh->cells_count ];

    vtkIdType cellId, faceId;
    int cellType;
    //vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
	for( cellId=0; cellId<mesh->cells_count; cellId++)
	{
        cellType =  vtkGrid->GetCellType(cellId);
        if( cellType == VTK_TETRA )
        {
          vtkTetra *cell = vtkTetra::SafeDownCast( vtkGrid->GetCell(cellId) );
          if( cell ) {
              mesh->cells[cellId].indices[0] = cell->GetPointId(0);
              mesh->cells[cellId].indices[1] = cell->GetPointId(1);
              mesh->cells[cellId].indices[2] = cell->GetPointId(2);
              mesh->cells[cellId].indices[3] = cell->GetPointId(3);
            }
        }
        else if( cellType == VTK_QUADRATIC_TETRA )
        {
          vtkQuadraticTetra *cell = vtkQuadraticTetra::SafeDownCast( vtkGrid->GetCell(cellId) );
          if( cell ) {
              mesh->cells[cellId].indices[0] = cell->GetPointId(0);
              mesh->cells[cellId].indices[1] = cell->GetPointId(1);
              mesh->cells[cellId].indices[2] = cell->GetPointId(2);
              mesh->cells[cellId].indices[3] = cell->GetPointId(3);
            }
        } else {
          // copy point ids
          for (int j = 0; j < 4; j++)
              mesh->cells[cellId].indices[j] = -1;
        }
	}

	printf("done.\n");

	printf("[VTKLoader] Loading mesh point attributes..");

	vtkDataArray* da = data->GetVectors();

    if( !da )
        da = data->GetArray("velocity");
	if( !da )
		da = data->GetArray("momentum");

	if( !da || da->GetNumberOfComponents() < 3)
	{
		("\nPoint data could not be loaded! Available data: ");
		int numArrays = data->GetNumberOfArrays();
		for( int i=0; i<numArrays; i++){
			const char* name = data->GetArrayName(i);
			printf("%s;",name);
		}
	}
	else
	{
		// Copy Point data
		delete mesh->node_attributes;
		mesh->node_attributes = new VertexAttribute[ mesh->node_count ];

		for( int i=0; i< mesh->node_count; i++)
		{
			double *v = da->GetTuple3(i);

			mesh->node_attributes[i].x = (float)v[0] * 100.0f;
			mesh->node_attributes[i].y = (float)v[1] * 100.0f;
			mesh->node_attributes[i].z = (float)v[2] * 100.0f;
			mesh->node_attributes[i].w = (float) 0.0 * 100.0f;
		}
	}

	printf("done.\n");
}
