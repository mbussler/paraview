#include "vtkTetRidge2D.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLine.h"
#include "vtkMergePoints.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTriangle.h"
#include "vtkIncrementalPointLocator.h"
#include "vtkTransform.h"

#include <math.h>
#include <algorithm>
#include "vtkUnstructuredGrid.h"

#include "linalg.h"
#include <map>
#include <set>

vtkStandardNewMacro(vtkTetRidge2D);

inline void swap( int& a, int& b) {
    int swp=a; a=b; b=swp;
};

//----------------------------------------------------------------------------
vtkTetRidge2D::vtkTetRidge2D()
{
    // by default process active point scalars
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::SCALARS);

    // by default process active point vectors
    this->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::VECTORS);

    // by default process active point vectors
    this->SetInputArrayToProcess(2,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::VECTORS);

    this->IsoValue = 0.5;
}

//----------------------------------------------------------------------------
vtkTetRidge2D::~vtkTetRidge2D()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int vtkTetRidge2D::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
// get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkPolyData *input = vtkPolyData::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkIdType cellId, i, updateTime;
    vtkSmartPointer<vtkCellArray> newLines;
    vtkSmartPointer<vtkPoints> newPoints;

    vtkIdType estimatedSize, numCells=input->GetNumberOfCells();
    vtkIdType numPts=input->GetNumberOfPoints();
    vtkPoints *inPts=input->GetPoints();
    int numberOfPoints;
    vtkPointData *inPD=input->GetPointData(), *outPD = output->GetPointData();
    vtkCellData *inCD=input->GetCellData(), *outCD = output->GetCellData();
    vtkCellData *outClippedCD = NULL;

    vtkDataArray *scalarArray = this->GetInputArrayToProcess(0, inputVector);
    if( !scalarArray )
    {
        vtkErrorMacro("No input array selected");
        return 0;
    }
    vtkDataArray *gradients = this->GetInputArrayToProcess(1,inputVector);
    if( !gradients )
    {
        vtkErrorMacro("No input gradients selected");
        return 0;
    }
    vtkDataArray *eV = this->GetInputArrayToProcess(2,inputVector);
    if( !eV )
    {
        vtkErrorMacro("No input eigen vectors selected");
        return 0;
    }

    vtkDebugMacro(<< "Calculating iso surface");

    // Initialize self; create output objects
    //
    if ( !input || numPts < 1 || inPts == NULL )
    {
        vtkDebugMacro(<<"No data.");
        return 1;
    }

    // new point data
    newPoints = vtkSmartPointer<vtkPoints>::New();
    newLines  = vtkSmartPointer<vtkCellArray>::New();

    // copy point data arrays
    int type = scalarArray->GetDataType();
    vtkSmartPointer<vtkDataArray> values = vtkDataArray::CreateDataArray(type);
    values->SetNumberOfComponents(scalarArray->GetNumberOfComponents());
    values->SetName(scalarArray->GetName());
    outPD->AddArray(values);

    // iterate over all cells, insert triangles forming the iso surface
    for( cellId=0; cellId<numCells; cellId++) 
    {
        vtkCell* cell = input->GetCell(cellId);
        if( cell && cell->GetCellType() == VTK_TRIANGLE ) 
        {
            vtkSmartPointer<vtkIdList> pointIds = cell->GetPointIds();

            int ptIds[3];
            ptIds[0] = pointIds->GetId(0);
            ptIds[1] = pointIds->GetId(1);
            ptIds[2] = pointIds->GetId(2);

            double ev[3*3], lev[3];
            eV->GetTuple(ptIds[0], &ev[0*3]);
            eV->GetTuple(ptIds[1], &ev[1*3]);
            eV->GetTuple(ptIds[2], &ev[2*3]);
            getLargestEV(ev, lev);

            double gradDotEV[3];
            for (int n=0; n<3; ++n) 
            {
                double orientedEV[3];
                orientedEV[0] = ev[n*3+0];
                orientedEV[1] = ev[n*3+1];
                orientedEV[2] = ev[n*3+2];

                double grads[3];
                gradients->GetTuple(ptIds[n], grads);

                if( this->Valley ){
                    vec3scal(grads, -1.0, grads);
                }

                if (vec3dot( lev, orientedEV) < 0) {
                    vec3scal(orientedEV, -1.0, orientedEV);
                }
                gradDotEV[n] = vec3dot(grads, orientedEV);
            }

            int lineindex = 0;

            XYZ p0, p1, p2;
            input->GetPoint(ptIds[0], &p0.x);
            input->GetPoint(ptIds[1], &p1.x);
            input->GetPoint(ptIds[2], &p2.x);

            double v0, v1, v2;
            scalarArray->GetTuple(ptIds[0], &v0);
            scalarArray->GetTuple(ptIds[1], &v1);
            scalarArray->GetTuple(ptIds[2], &v2);

            if( gradDotEV[0] < 0 ) lineindex |= 1;
            if( gradDotEV[1] < 0 ) lineindex |= 2;
            if( gradDotEV[2] < 0 ) lineindex |= 4;

            XYZ np0, np1;
            double val0, val1;

            switch( lineindex )
            {
            case 0x00:
            case 0x07:
                break;
            case 0x06:
            case 0x01:                
                np0 = VertexInterpolate(0.0,p0,p1,gradDotEV[0],gradDotEV[1], v0, v1, val0);
                np1 = VertexInterpolate(0.0,p0,p2,gradDotEV[0],gradDotEV[2], v0, v2, val1);
                if((val0 > MinDataValue && val0 < MaxDataValue ) &&
                   (val1 > MinDataValue && val1 < MaxDataValue ) ) 
                {
                    InsertLine(np0, np1, newPoints, newLines);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                }
                break;
            case 0x05:
            case 0x02:
                np0 = VertexInterpolate(0.0,p1,p0,gradDotEV[1],gradDotEV[0], v1, v0, val0);
                np1 = VertexInterpolate(0.0,p1,p2,gradDotEV[1],gradDotEV[2], v1, v2, val1);
                if((val0 > MinDataValue && val0 < MaxDataValue ) &&
                   (val1 > MinDataValue && val1 < MaxDataValue ) )
                {
                    InsertLine(np0, np1, newPoints, newLines);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                }
                break;
            case 0x04:
            case 0x03:
                np0 = VertexInterpolate(0.0,p2,p0,gradDotEV[2],gradDotEV[0], v2, v0, val0);
                np1 = VertexInterpolate(0.0,p2,p1,gradDotEV[2],gradDotEV[1], v2, v1, val1);
                if((val0 > MinDataValue && val0 < MaxDataValue ) &&
                   (val1 > MinDataValue && val1 < MaxDataValue ) )
                {
                    InsertLine(np0, np1, newPoints, newLines);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                }
                break;
            }
        }

        this->UpdateProgress((float)cellId/numCells);
    }

    //std::vector<unsigned int> indices_merged;
    //std::vector<XYZ> vertices_merged;
    //std::vector<double> scalars_merged;
    //mergeTriangles(newPoints, values, indices_merged, vertices_merged, scalars_merged);

    //std::vector<int> triLabels;
    //std::vector<int> compSizes;
    //labelComponents(&indices_merged, Manifold, &triLabels, &compSizes);

    //std::cout << "Region triangle count threshold: " << RegionThreshold << std::endl;
    //std::set<int> drawIds;

    //int largestRegId = -1;
    //for (int nComp = 0; nComp < compSizes.size(); ++nComp) {
    //    if (compSizes[nComp] > RegionThreshold)
    //        drawIds.insert(nComp);
    //}

    //// std::cout << "Component #" << largestRegId << " has " << compSizes[largestRegId] << " triangles." << std::endl;
    //// std::cout << "triLabels.size(): " << triLabels.size() << std::endl;
    //std::cout << "indices_merged.size(): " << indices_merged.size() << std::endl;
    //std::cout << "vertices_merged.size(): " << vertices_merged.size() << std::endl;
    //std::cout << "scalars_merged.size(): " << scalars_merged.size() << std::endl;

    //newPoints =  vtkPoints::New();
    //newLines = vtkCellArray::New();
    //values->SetNumberOfTuples(0);

    //for (int t=0; t<triLabels.size(); ++t) {
    //    vtkIdType id[1];
    //    int vIdx; 
    //    vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();
    //    bool drawTri = (drawIds.find(triLabels[t]) != drawIds.end());
    //    if (drawTri) {

    //        vIdx = indices_merged[t*3+0];
    //        id[0] = newPoints->InsertNextPoint(&vertices_merged[vIdx].x);
    //        values->InsertNextTuple1(scalars_merged[vIdx]);
    //        tri->GetPointIds()->SetId(0, id[0]);

    //        vIdx = indices_merged[t*3+1];
    //        id[0] = newPoints->InsertNextPoint(&vertices_merged[vIdx].x);
    //        values->InsertNextTuple1(scalars_merged[vIdx]);
    //        tri->GetPointIds()->SetId(1, id[0]);

    //        vIdx = indices_merged[t*3+2];
    //        id[0] = newPoints->InsertNextPoint(&vertices_merged[vIdx].x);
    //        values->InsertNextTuple1(scalars_merged[vIdx]);
    //        tri->GetPointIds()->SetId(2, id[0]);

    //        newLines->InsertNextCell(tri);
    //    }
    //}

    output->SetPoints(newPoints);
    output->SetLines(newLines);

    output->Squeeze();

    return 1;
}



XYZ vtkTetRidge2D::VertexInterpolate( double isovalue, XYZ p1, XYZ p2, double valp1, double valp2, double data1, double data2, double& val)
{
    XYZ p;
    double mu;
    
    if( abs(isovalue-valp1) < 0.00001)
        return p1;
    if( abs(isovalue-valp2) < 0.00001)
        return p2;
    if( abs(valp1-valp2) < 0.00001)
        return p1;
    mu = (isovalue - valp1) / (valp2-valp1);
    p.x = p1.x + mu * (p2.x - p1.x);
    p.y = p1.y + mu * (p2.y - p1.y);
    p.z = p1.z + mu * (p2.z - p1.z);
    val = data1 + mu * (data2 - data1);

    return p;
}

void vtkTetRidge2D::InsertLine(XYZ np0, XYZ np1, vtkPoints * newPoints, vtkCellArray * newLines)
{
    newLines->InsertNextCell(2);
    newLines->InsertCellPoint( newPoints->InsertNextPoint(&np0.x));
    newLines->InsertCellPoint( newPoints->InsertNextPoint(&np1.x));
}

int vtkTetRidge2D::FillInputPortInformation(int port, vtkInformation* info)
{
    if ( port == 0 ) {
        info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData" );
        return 1;
    }
    return 0;
}


void vtkTetRidge2D::getLargestEV(double *ev, double lev[3])
{
    // calculate largest eigen vector of cell
    double evSet[3][3*2];
    for (int n=0; n<3; ++n) {
        evSet[0][2*n+0] =  ev[3*n+0];
        evSet[1][2*n+0] =  ev[3*n+1];
        evSet[2][2*n+0] =  ev[3*n+2];
        evSet[0][2*n+1] = -ev[3*n+0];
        evSet[1][2*n+1] = -ev[3*n+1];
        evSet[2][2*n+1] = -ev[3*n+2];
    }
    //Covariance Matrix
    mat3 C;
    for(int j=0; j<3; ++j) {
        for(int i=0; i<3; ++i) {
            double sum = 0.0;
            for(int n=0; n<3; ++n){
                sum += evSet[i][n] * evSet[j][n];
            }
            C[j][i] = sum/9.0;
        }
    }
    vec3 lambda;
    mat3eigenvalues( C, lambda);

    // get EV of largest EW
    int largest = 0;
    if( lambda[0] < lambda[1]) largest = 1;
    if( lambda[1] < lambda[2]) largest = 2;
    mat3realEigenvector(C, lambda[largest], lev);
}

void vtkTetRidge2D::mergeTriangles( vtkPoints* points, vtkDataArray* values, 
    std::vector<unsigned int>& indices, 
    std::vector<XYZ>& vertices, 
    std::vector<double>& scalars )
{

    indices.resize(0);
    vertices.resize(0);

    std::map<FKeys3, unsigned int> vertexInfo;

    unsigned int vertexID = 0;

    // go over all triangles
    for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ptId++) {

        // generate vertex key
        XYZ pos;
        points->GetPoint(ptId, &pos.x);
        FKeys3 key(pos.x, pos.y, pos.z);

        // find vertex
        if (vertexInfo.find(key) == vertexInfo.end()) {
            // not found

            vertexInfo[key] = vertexID;

            // insert vertex and its corresponding data
            vertices.push_back(pos);
            scalars.push_back(values->GetTuple1(ptId));

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

void vtkTetRidge2D::labelComponents(std::vector<unsigned int>* indices, 
    bool manifold,
    std::vector<int> *triangleLabels,
    std::vector<int> *objectSizes)
{ 
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

                compSize += labelComponents_rek( indices, 
                                                 manifold, 
                                                 t, 
                                                 &edges, 
                                                 triaVisited, 
                                                 &flippedCnt, 
                                                 connCompCnt, 
                                                 triangleLabels);

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

int vtkTetRidge2D::labelComponents_rek(std::vector<unsigned int>* indices,
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