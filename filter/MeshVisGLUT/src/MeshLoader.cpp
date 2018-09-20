#include "MeshLoader.h"

// file I/O
#include <iostream>
#include <fstream>
using namespace std;

MeshLoader::MeshLoader()
{
};

MeshLoader::~MeshLoader()
{
};

bool MeshLoader::loadNeighborsFromFile( const char* filename, TetMeshPtr mesh)
{
	char neighfn[255];
	sprintf(neighfn,"%s.nei",filename);

	printf("[MeshLoader] Loading neighbors from %s..", neighfn);

	ifstream neighFile ( neighfn, ios::binary );
	if( neighFile.is_open() )
	{
		neighFile.seekg(0);
		neighFile.read((char*)&mesh->cells_count, sizeof(int));

		mesh->neighbors = new Cell[ mesh->cells_count ];

		list<Face> outer_faces;
		list<Vector> normals;
		list<Vector> face_vertices;

		/* read binary */
		for( int i=0; i<mesh->cells_count; i++ )
		{
			for(int j=0; j<4; j++)
			{
				neighFile.read((char*)&mesh->neighbors[i].indices[j], sizeof(int));

				if( mesh->neighbors[i].indices[j] == -1 )
				{
					Face f;
					f.pts[0] = mesh->cells[i].indices[ faces[j][0] ];
					f.pts[1] = mesh->cells[i].indices[ faces[j][1] ];
					f.pts[2] = mesh->cells[i].indices[ faces[j][2] ];

					// calculate face normal
					Vector p1 = Pt2Vec( mesh->nodes[f.pts[0]] );
					Vector p2 = Pt2Vec( mesh->nodes[f.pts[1]] );
					Vector p3 = Pt2Vec( mesh->nodes[f.pts[2]] );

					Vector f_center = 1.0f/3.0f * ( p1 + p2 + p3);
					Vector normal = VectorNormalize(VectorCross( p1-p2, p2-p3));

					// calculate cell center
					Cell c = mesh->cells[i];
					Vector vs[4];

					for( int j=0; j<4; j++) {
						vs[j] = Pt2Vec( mesh->nodes[ mesh->cells[i].indices[j] ] );
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
			}
		}
		neighFile.close();

		mesh->copyOuterFaces(outer_faces);
		mesh->copyNormals(normals);
		mesh->copyFaceVertices(face_vertices);

		printf("done.\n");
		return true;
	}
	else
	{
		printf("file opening failed!\n");
		return false;
	}
};

void MeshLoader::saveNeighborsToFile( const char* filename, TetMeshPtr mesh )
{
	char neighfn[255];
	sprintf(neighfn,"%s.nei",filename);

	printf("[MeshLoader] Saving neighbors to %s..", neighfn);

	ofstream neighFile ( neighfn, ios::binary | ios::trunc );

	if( neighFile.is_open() )
	{
		neighFile.seekp(0);
		neighFile.write((char*)&mesh->cells_count, sizeof(int));

		/* write binary */
		for(int i=0; i<mesh->cells_count; i++)
		{
			for(int j=0; j<4; j++)
			{
				neighFile.write((char*)&mesh->neighbors[i].indices[j], sizeof(int));
			}
		}

		neighFile.close();
		printf("done.\n");
	}
	else
	{
		printf("error: could not open %s for write.\n", neighfn);
	}
};
