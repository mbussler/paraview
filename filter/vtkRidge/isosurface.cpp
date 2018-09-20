#include "isosurface.h"
#include <set>
#include <map>

//---Isosurface---------------------------------------------------------------
Isosurface::Isosurface() {
    m_VTKdataPoints = vtkSmartPointer<vtkPoints>::New();    
    m_VTKisoTriangles = vtkSmartPointer<vtkCellArray>::New();    
    m_VTKPolyData = vtkSmartPointer<vtkPolyData>::New();
    m_VTKdataGrads = vtkSmartPointer<vtkFloatArray>::New();
    m_VTKisoGrads = vtkSmartPointer<vtkFloatArray>::New();
    m_VTKPointData = vtkSmartPointer<vtkDoubleArray>::New();
    m_VTKPointData->SetNumberOfComponents(1);
    m_VTKPointData->SetName("Density");
    m_EigenValues = vtkSmartPointer<vtkImageData>::New();
    m_dim[0] = 0;
    m_dim[1] = 0;
    m_dim[2] = 0;
    m_data = NULL;
    m_grad_data = NULL;
    m_hesse_data = NULL;
    m_hasGradients = false;
    m_origin[0] = 0.0;
    m_origin[1] = 0.0;
    m_origin[2] = 0.0;
    m_valley = false;
}

Isosurface::~Isosurface(){
    CleanUp();
}

void Isosurface::generateTriangles(int regThreshold, double angleThreshold, double minDataValue, double maxDataValue, double featureThreshold, double waveSpeed, unsigned int stencilRange){

    // double *vertsCPU = NULL;
    double *verts = NULL;
    double *evals = NULL;
    double *pointData = NULL;
    std::set<unsigned int> nodeIndicesTmp;

    uint numVertsGPU;

    std::cout << featureThreshold << std::endl;
    if (!generateRidgeSurface(&numVertsGPU, &verts, &evals, &pointData, featureThreshold, stencilRange))
    {
        return;
    }

    cout << numVertsGPU << " vertices before filtering." << endl;

    std::vector< float3 > angleFilteredSoup;
    std::vector< float > angleFilteredData;

    for(int i=0; i<numVertsGPU; i+=3) 
    {
        vtkIdType id[1];
        vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();

        //double* vertexCoord = &verts[(i+0)*3];

        //vec3 vertex0to1 = {vertexCoord[3] - vertexCoord[0],
        //    vertexCoord[4] - vertexCoord[1],
        //    vertexCoord[5] - vertexCoord[2]};

        //vec3 vertex0to2 = {vertexCoord[6] - vertexCoord[0],
        //    vertexCoord[7] - vertexCoord[1],
        //    vertexCoord[8] - vertexCoord[2]};

        //vec3 normal;
        //vec3cross(vertex0to1, vertex0to2, normal);
        //vec3nrm(normal,normal);

        //vec3 filterDirection = {0.0, 0.0, 1.0};
        //vec3nrm(filterDirection, filterDirection);

        //if (vec3dot(normal, filterDirection) < 0)
        //{
        //    vec3scal(normal, -1, normal);
        //}

        ////double waveSpeed = 2.0;
        //double cosWaveAngle = waveSpeed / sqrt(1.0 + waveSpeed*waveSpeed);
        //double waveAngle = acos(cosWaveAngle) * 180 / M_PI;

        //double normalToWaveAngle = acos(vec3dot(normal,filterDirection)) * 180 / M_PI;
        ////cout << waveAngle << " " << normalToWaveAngle << " " << abs(waveAngle - normalToWaveAngle) << endl;
        //if (abs(waveAngle - normalToWaveAngle) > angleThreshold){
        //    continue;
        //}

        if (minDataValue < maxDataValue){
            if (pointData[i] < minDataValue || pointData[i] > maxDataValue) continue;
            if (pointData[i+1] < minDataValue || pointData[i+1] > maxDataValue) continue;
            if (pointData[i+2] < minDataValue || pointData[i+2] > maxDataValue) continue;
        } else {
            if (pointData[i] < minDataValue && pointData[i] > maxDataValue) continue;
            if (pointData[i+1] < minDataValue && pointData[i+1] > maxDataValue) continue;
            if (pointData[i+2] < minDataValue && pointData[i+2] > maxDataValue) continue;
        }

        if (regThreshold) {
            float3 v;
            v.x = verts[(i+0)*3+0]; v.y = verts[(i+0)*3+1]; v.z = verts[(i+0)*3+2];
            angleFilteredSoup.push_back(v);
            angleFilteredData.push_back(pointData[i]);
            v.x = verts[(i+1)*3+0]; v.y = verts[(i+1)*3+1]; v.z = verts[(i+1)*3+2];
            angleFilteredSoup.push_back(v);
            angleFilteredData.push_back(pointData[i+1]);
            v.x = verts[(i+2)*3+0]; v.y = verts[(i+2)*3+1]; v.z = verts[(i+2)*3+2];
            angleFilteredSoup.push_back(v);
            angleFilteredData.push_back(pointData[i+2]);
        }
        else {
            id[0] = m_VTKdataPoints->InsertNextPoint(&verts[(i+0)*3]);
            // m_VTKisoGrads->InsertTuple(id[0], Isosurface::normalize(&trigrads[(i+0)*3]));
            tri->GetPointIds()->SetId(0, id[0]);

            id[0] = m_VTKdataPoints->InsertNextPoint(&verts[(i+1)*3]);
            // m_VTKisoGrads->InsertTuple(id[0], Isosurface::normalize(&trigrads[(i+1)*3]));
            tri->GetPointIds()->SetId(1, id[0]);

            id[0] = m_VTKdataPoints->InsertNextPoint(&verts[(i+2)*3]);
            // m_VTKisoGrads->InsertTuple(id[0], Isosurface::normalize(&trigrads[(i+2)*3]));
            tri->GetPointIds()->SetId(2, id[0]);

            m_VTKisoTriangles->InsertNextCell(tri);

            m_VTKPointData->InsertNextValue(pointData[i]);
            m_VTKPointData->InsertNextValue(pointData[i+1]);
            m_VTKPointData->InsertNextValue(pointData[i+2]);
        }

        // size_t nodeIndex = 0;
        // nodeIndex += static_cast<size_t>(floor((vertexCoord[0] + vertexCoord[3] + vertexCoord[6])/3.0 - m_origin[0])/m_spacing[0]);
        // nodeIndex += static_cast<size_t>(floor((vertexCoord[1] + vertexCoord[4] + vertexCoord[7])/3.0 - m_origin[1])/m_spacing[1]) * m_dim[0];
        // nodeIndex += static_cast<size_t>(floor((vertexCoord[2] + vertexCoord[5] + vertexCoord[8])/3.0 - m_origin[2])/m_spacing[2]) * m_dim[0] * m_dim[1];

        // nodeIndicesTmp.insert(nodeIndex);
        // nodeIndicesTmp.insert(nodeIndex + 1);
        // nodeIndicesTmp.insert(nodeIndex + m_dim[0]);
        // nodeIndicesTmp.insert(nodeIndex + m_dim[0] + 1);
        // nodeIndicesTmp.insert(nodeIndex + m_dim[1]*m_dim[0]);
        // nodeIndicesTmp.insert(nodeIndex + m_dim[1]*m_dim[0] + 1);
        // nodeIndicesTmp.insert(nodeIndex + m_dim[1]*m_dim[0] + m_dim[0]);
        // nodeIndicesTmp.insert(nodeIndex + m_dim[1]*m_dim[0] + m_dim[0] + 1);
    }

    ////////////////////////////////////////////THE SOUP IS SERVED
    if (regThreshold)
    {
        std::cout << "Angle filtered soup contains: " << angleFilteredSoup.size() << " vertices" << std::endl;

        std::vector<unsigned int> indices_merged;
        std::vector<float3> vertices_merged;
        std::vector<float> scalars_merged;
        std::vector<int> triLabels;
        std::vector<int> compSizes;

        indices_merged.clear();
        vertices_merged.clear();
        scalars_merged.clear();

        mergeTriangles(angleFilteredSoup.data(), angleFilteredData.data(), angleFilteredSoup.size(), indices_merged, vertices_merged, scalars_merged);
        labelComponents(&indices_merged, m_manifold, &triLabels, &compSizes);

        std::cout << "Region triangle count threshold: " << regThreshold << std::endl;
        std::set<int> drawIds;

        int largestRegId = -1;
        for (int nComp = 0; nComp < compSizes.size(); ++nComp) {
            if (compSizes[nComp] > regThreshold)
                drawIds.insert(nComp);
        }

        // std::cout << "Component #" << largestRegId << " has " << compSizes[largestRegId] << " triangles." << std::endl;
        // std::cout << "triLabels.size(): " << triLabels.size() << std::endl;
        std::cout << "indices_merged.size(): " << indices_merged.size() << std::endl;
        std::cout << "vertices_merged.size(): " << vertices_merged.size() << std::endl;
        std::cout << "scalars_merged.size(): " << scalars_merged.size() << std::endl;

        int test = 0;
        for (int t=0; t<triLabels.size(); ++t) {
            vtkIdType id[1];
            int vIdx; 
            vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();
            bool drawTri = (drawIds.find(triLabels[t]) != drawIds.end());
            if (drawTri) {
                ++test;

                vIdx = indices_merged[t*3+0];
                id[0] = m_VTKdataPoints->InsertNextPoint((float*)&vertices_merged[vIdx]);
                m_VTKPointData->InsertNextValue(scalars_merged[vIdx]);
                tri->GetPointIds()->SetId(0, id[0]);

                vIdx = indices_merged[t*3+1];
                id[0] = m_VTKdataPoints->InsertNextPoint((float*)&vertices_merged[vIdx]);
                m_VTKPointData->InsertNextValue(scalars_merged[vIdx]);
                tri->GetPointIds()->SetId(1, id[0]);

                vIdx = indices_merged[t*3+2];
                id[0] = m_VTKdataPoints->InsertNextPoint((float*)&vertices_merged[vIdx]);
                m_VTKPointData->InsertNextValue(scalars_merged[vIdx]);
                tri->GetPointIds()->SetId(2, id[0]);

                m_VTKisoTriangles->InsertNextCell(tri);
            }
        }
    }
    ////////////////////////////////////////////BON APPETIT!

    m_EigenValues->SetDimensions(m_dim[0], m_dim[1], m_dim[2]);
    m_EigenValues->SetSpacing(m_spacing);
    m_EigenValues->SetExtent(m_extent);
    m_EigenValues->SetOrigin(m_origin);
    m_EigenValues->AllocateScalars(VTK_DOUBLE, 3);
    double *imageData = static_cast<double *>(m_EigenValues->GetScalarPointer());
    memcpy(imageData, evals, m_dim[0]*m_dim[1]*m_dim[2]*3*sizeof(double));

    if (verts) delete[] verts; verts = NULL;
    if (evals) delete[] evals; evals = NULL;
    if (pointData) delete [] pointData; pointData = NULL;
    CleanUp();

    //    for (size_t i = 0; i < numNodeIndices; i++)
    //    {
    //        vtkIdType id[1];
    //        float x, y, z;
    //        x = nodeIndices[i] % m_dim[0];
    //        y = (nodeIndices[i] / m_dim[0]) % m_dim[1];
    //        z = nodeIndices[i] / (m_dim[0] * m_dim[1]);

    //        id[0] = m_VTKdataPoints->InsertNextPoint(x, y, z);;
    //        vertices->InsertNextCell(1, id);
    //    }

    //    // TODO: account for scale!
    //    for (std::set<unsigned int>::iterator it = nodeIndicesTmp.begin(); it != nodeIndicesTmp.end(); it++)
    //    {
    //        vtkIdType id[1];
    //        float x, y, z;
    //        x = m_spacing[0]*(*it % m_dim[0]) + m_origin[0];
    //        y = m_spacing[1]*((*it / m_dim[0]) % m_dim[1]) + m_origin[1];
    //        z = m_spacing[2]*(*it / (m_dim[0] * m_dim[1])) + m_origin[2];

    //        id[0] = m_RidgeNodesDataPoints->InsertNextPoint(x, y, z);;
    //        m_RidgeNodesVertices->InsertNextCell(1, id);
    //        m_RidgeNodesIndices->InsertNextValue(*it);
    //    }

    // m_VTKPolyData->GetPointData()->SetNormals(m_VTKisoGrads);
    cout << "Storing " << m_VTKdataPoints->GetNumberOfPoints() << " points." << endl;
    m_VTKPolyData->SetPoints(m_VTKdataPoints);

    cout << "Storing " << m_VTKisoTriangles->GetNumberOfCells() << " triangles." << endl;
    m_VTKPolyData->SetPolys(m_VTKisoTriangles);

    // vtkSmartPointer<vtkFloatArray> pointDataVtk = vtkSmartPointer<vtkFloatArray>::New();
    // vtkSmartPointer<vtkDoubleArray> pointDataVtk = vtkSmartPointer<vtkDoubleArray>::New();
    // pointDataVtk->SetArray(m_pointDataRaw.data(), m_VTKdataPoints->GetNumberOfPoints(), 1);
    // pointDataVtk->SetArray(pointData, m_VTKdataPoints->GetNumberOfPoints(), 1);
    cout << "Storing " << m_VTKPointData->GetNumberOfTuples() << " values." << endl;
    m_VTKPolyData->GetPointData()->AddArray(m_VTKPointData);

    //m_VTKPolyData->GetPointData()->AddArray(evalsArray);
    //evalsArray->Delete();

    //    m_RidgeNodesPolyData->SetPoints(m_RidgeNodesDataPoints);
    //    m_RidgeNodesPolyData->SetVerts(m_RidgeNodesVertices);
    //    m_RidgeNodesPolyData->GetPointData()->AddArray(m_RidgeNodesIndices);


    cout << "Triangles generated. Exiting method." << endl;
}

//----Ridge surface-----------------------------------------------------------------------------

bool Isosurface::generateRidgeSurface(uint* numVerts, double** verts, double **evals, double **pointData, double featureThreshold, unsigned int stencilRange) {
    cout << "\nRidge surface on GPU" << endl;
    double ridgeOrigin[3] = {m_origin[0] + m_spacing[0] * static_cast<double>(m_extent[0]),
        m_origin[1] + m_spacing[1] * static_cast<double>(m_extent[2]),
        m_origin[2] + m_spacing[2] * static_cast<double>(m_extent[4])};
    return cudaRidgeData( m_dim, m_data, m_grad_data, m_hesse_data,
        evals, numVerts, verts, pointData, 
        featureThreshold, ridgeOrigin, m_spacing, stencilRange, m_valley);
}


//----Cleanup-----------------------------------------------------------------------------------
void Isosurface::CleanUp(void) {
    m_data=NULL;
}


void Isosurface::GetEdgeIndices(int index, int* outArray ){
    outArray[0*2+0] = index;
    outArray[0*2+1] = index+1;

    outArray[1*2+0] = index+1;
    outArray[1*2+1] = index+1+m_dim[0];

    outArray[2*2+0] = index+m_dim[0];
    outArray[2*2+1] = index+1+m_dim[0];

    outArray[3*2+0] = index;
    outArray[3*2+1] = index+m_dim[0];

    outArray[4*2+0] = index+m_dim[0]*m_dim[1];
    outArray[4*2+1] = index+1+m_dim[0]*m_dim[1];

    outArray[5*2+0] = index+1+m_dim[0]*m_dim[1];
    outArray[5*2+1] = index+1+m_dim[0]+m_dim[0]*m_dim[1];

    outArray[6*2+0] = index+m_dim[0]+m_dim[0]*m_dim[1];
    outArray[6*2+1] = index+1+m_dim[0]+m_dim[0]*m_dim[1];

    outArray[7*2+0] = index+m_dim[0]*m_dim[1];
    outArray[7*2+1] = index+m_dim[0]+m_dim[0]*m_dim[1];

    outArray[8*2+0] = index;
    outArray[8*2+1] = index+m_dim[0]*m_dim[1];

    outArray[9*2+0] = index+1;
    outArray[9*2+1] = index+1+m_dim[0]*m_dim[1];

    outArray[10*2+0] = index+1+m_dim[0];
    outArray[10*2+1] = index+1+m_dim[0]+m_dim[0]*m_dim[1];

    outArray[11*2+0] = index+m_dim[0];
    outArray[11*2+1] = index+m_dim[0]+m_dim[0]*m_dim[1];

}

void Isosurface::GetEdgeVertices(double* verts, double* pVerts[12] ) {
    verts[0*3+0]=double(0.5); pVerts[0]=&verts[0*3+0];
    verts[0*3+1]=double(0);
    verts[0*3+2]=double(0);

    verts[1*3+0]=double(1);
    verts[1*3+1]=double(0);
    verts[1*3+2]=double(0.5); pVerts[1]=&verts[1*3+2];

    verts[2*3+0]=double(0.5); pVerts[2]=&verts[2*3+0];
    verts[2*3+1]=double(0);
    verts[2*3+2]=double(1); 

    verts[3*3+0]=double(0);
    verts[3*3+1]=double(0);
    verts[3*3+2]=double(0.5); pVerts[3]=&verts[3*3+2];

    verts[4*3+0]=double(0.5); pVerts[4]=&verts[4*3+0];
    verts[4*3+1]=double(1);
    verts[4*3+2]=double(0);

    verts[5*3+0]=double(1);
    verts[5*3+1]=double(1);
    verts[5*3+2]=double(0.5); pVerts[5]=&verts[5*3+2];

    verts[6*3+0]=double(0.5); pVerts[6]=&verts[6*3+0];
    verts[6*3+1]=double(1);
    verts[6*3+2]=double(1);

    verts[7*3+0]=double(0);
    verts[7*3+1]=double(1);
    verts[7*3+2]=double(0.5); pVerts[7]=&verts[7*3+2];

    verts[8*3+0]=double(0);
    verts[8*3+1]=double(0.5); pVerts[8]=&verts[8*3+1];
    verts[8*3+2]=double(0);

    verts[9*3+0]=double(1);
    verts[9*3+1]=double(0.5); pVerts[9]=&verts[9*3+1];
    verts[9*3+2]=double(0);

    verts[10*3+0]=double(1);
    verts[10*3+1]=double(0.5); pVerts[10]=&verts[10*3+1];
    verts[10*3+2]=double(1);

    verts[11*3+0]=double(0);
    verts[11*3+1]=double(0.5); pVerts[11]=&verts[11*3+1];
    verts[11*3+2]=double(1);
}

void  Isosurface::mergeTriangles(float3* pos,
    float* data,
    int totalVerts, 
    std::vector<unsigned int>& indices, 
    std::vector<float3>& vertices,
    std::vector<float>& scalars)
{ // indices 3 indices per triangle

    indices.resize(0);
    vertices.resize(0);
    // normals.resize(0);
    // samples.resize(0);

    std::map<FKeys3, unsigned int> vertexInfo;

    unsigned int vertexID = 0;

    // go over all triangles
    for (int t = 0; t < totalVerts/3; t++) {

        // for each vertex in triangle
        for (int v = 0; v < 3; v++) {

            // generate vertex key
            FKeys3 key(pos[t*3+v].x, pos[t*3+v].y, pos[t*3+v].z);

            // find vertex
            if (vertexInfo.find(key) == vertexInfo.end()) {
                // not found

                vertexInfo[key] = vertexID;

                // insert vertex and its corresponding data
                vertices.push_back(pos[t*3+v]);
                scalars.push_back(data[t*3+v]);

                // insert indices
                indices.push_back(vertexID);
                vertexID++;
            }
            else {

                // found

                // insert indices
                indices.push_back(vertexInfo[key]);	
            }
        }
    }
}

int  Isosurface::labelComponents_rek(std::vector<unsigned int>* indices,
    bool manifold,
    int triangle,
    std::map<std::pair<int, int>, std::vector<int> > *edges,
    bool *triaVisited, int *flippedCnt,
    int compLabel,
    std::vector<int> *triangleComponents)
{
    int visitedCnt = 0;
    int triaCnt = indices->size() / 3;

    std::vector<int> stack;
    stack.push_back(triangle);

    while( !stack.empty() ) {

        triangle = stack.back();
        stack.pop_back();

        // go over edges of triangle
        for (int eIdx=0; eIdx<3; eIdx++) {

            int locIdx1 = eIdx;
            int locIdx2 = ( eIdx+1 < 3 ? eIdx+1 : 0 );

            // get vertex IDs
            int id1 = (*indices)[triangle*3 + locIdx1];
            int id2 = (*indices)[triangle*3 + locIdx2];

            // get lowId, highId
            int lowId = id1;
            int highId = id2;
            if (id1 > id2) {
                lowId = id2;
                highId = id1;
            }

            // get key
            std::pair<int, int> key;
            key.first = lowId;
            key.second = highId;

#if 0
            if ((*edges)[key].size() != 2) {
                // if ((*edges)[key].size() != 1) {
                //   printf("skipping edge shared by more than 2 triangles\n");
                // }
                continue;
            }
#else
            if ((*edges)[key].size() != 2) {
                if (manifold) 
                {
                    if ((*edges)[key].size() != 1) {
                        printf("skipping edge shared by more than 2 triangles %d\n", (*edges)[key].size());
                    }
                    continue;
                }
                else 
                {
                    if ((*edges)[key].size() == 1) 
                    {
                        continue;
                    }
                    else 
                    {
                        printf("edge shared by more than 2 triangles\n");
                    }
                }
            }
#endif

            // go over triangles of edge
            for (int tIdx=0; tIdx<(int)(*edges)[key].size(); tIdx++) {
                int nt = (*edges)[key][tIdx];
                if (!triaVisited[nt]) {

                    triaVisited[nt] = true;
                    visitedCnt++;
                    (*triangleComponents)[nt] = compLabel;

                    // push to stack
                    stack.push_back(nt);

                }
            }
        }


    }






    return visitedCnt;
}

void Isosurface::labelComponents(std::vector<unsigned int>* indices, 
    bool manifold,
    std::vector<int> *triangleLabels,
    std::vector<int> *objectSizes)
{ // triangleLabels, objectSizes: output
    int triaCnt = indices->size() / 3;

    triangleLabels->clear();
    triangleLabels->resize(triaCnt);
    objectSizes->clear();

    // get triangles that share each edge
    std::map<std::pair<int, int>, std::vector<int> > edges;  // edge, triangles
    for (int t=0; t<triaCnt; t++) {

        // go over edges of triangle
        for (int e=0; e<3; e++) {
            int locIdx1 = e;
            int locIdx2 = ( e+1 < 3 ? e+1 : 0 );

            // get vertex IDs
            int id1 = (*indices)[t*3 + locIdx1];
            int id2 = (*indices)[t*3 + locIdx2];

            // get lowId, highId
            int lowId = id1;
            int highId = id2;
            if (id1 > id2) {
                lowId = id2;
                highId = id1;
            }

            // collect info
            std::pair<int, int> key;
            key.first = lowId;
            key.second = highId;
            if (edges.find(key) == edges.end()) {
                // new edge
                std::vector<int> trias;
                trias.push_back(t);
                edges[key] = trias;
            }
            else {
                // existing edge
                edges[key].push_back(t);
            }
        }
        // go over vertices of triangle
        for (int e=0; e<3; e++) {
            int locIdx1 = e;

            // get vertex IDs
            int id1 = (*indices)[t*3 + locIdx1];

            // get lowId, highId
            int lowId = id1;
            int highId = id1;

            // collect info
            std::pair<int, int> key;
            key.first = lowId;
            key.second = highId;
            if (edges.find(key) == edges.end()) {
                // new edge
                std::vector<int> trias;
                trias.push_back(t);
                edges[key] = trias;
            }
            else {
                // existing edge
                edges[key].push_back(t);
            }
        }
    }

    //  printf("DEBUG> early return labelComponents\n");
    //  return ;

    // flip inconsistent triangles (if possible due to non-orientable surfaces)
    {  
        bool *triaVisited = new bool[triaCnt];
        for (int i=0; i<triaCnt; i++) { triaVisited[i] = false; };

        // loop over components
        int connCompCnt = 0;
        int flippedCnt = 0;
        //for (int t=0; t<triaCnt; t++) {
        for (int t=triaCnt-1; t>=0; t--) {

            if (!triaVisited[t]) {     
                // triangle defines orientation for current connected component
                triaVisited[t] = true;
                int compSize = 1;
                (*triangleLabels)[t] = connCompCnt;

                compSize += labelComponents_rek(indices, manifold, t, &edges, triaVisited, &flippedCnt, connCompCnt, triangleLabels);

                objectSizes->push_back(compSize);
                //printf("compSize = %d\n", compSize);
                connCompCnt++;
            }
        }

        // ####
        int visitedCnt = 0;
        for (int i=0; i<triaCnt; i++) {
            if (triaVisited[i]) {
                visitedCnt++;
            }
        }
        //    printf("%d visited and %d flipped triangles\n", visitedCnt, flippedCnt);
        delete[] triaVisited;
    }
}
