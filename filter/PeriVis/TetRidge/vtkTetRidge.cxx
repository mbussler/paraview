#include "vtkTetRidge.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
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

vtkStandardNewMacro(vtkTetRidge);

/* performance measure */
#include "timer.h"
#include <QElapsedTimer>

template <typename T>
inline void swap( T& a, T& b) {
    T swp=a; a=b; b=swp;
};

//----------------------------------------------------------------------------
vtkTetRidge::vtkTetRidge()
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

    // by default process active point scalars
    this->SetInputArrayToProcess(3,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::SCALARS);

    // by default process active point scalars
    this->SetInputArrayToProcess(4,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::TENSORS);

    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(2);
}

//----------------------------------------------------------------------------
vtkTetRidge::~vtkTetRidge()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int vtkTetRidge::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
// get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);

    // get the input and output
    vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkUnstructuredGrid *outgrid = vtkUnstructuredGrid::SafeDownCast(
        outInfo1->Get(vtkDataObject::DATA_OBJECT()));

    vtkIdType cellId, i, updateTime;
    vtkSmartPointer<vtkCellArray> newTris;
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

    vtkDataArray *eW = this->GetInputArrayToProcess(3,inputVector);
    if( !eW )
    {
        vtkErrorMacro("No input eigen values selected");
        return 0;
    }

    vtkDataArray *dataArray = this->GetInputArrayToProcess(4,inputVector);
    if( !dataArray )
    {
        vtkErrorMacro("No additional data values selected");
        return 0;
    }
    int numDataComponents = dataArray->GetNumberOfComponents();

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
    newTris   = vtkSmartPointer<vtkCellArray>::New();

    // copy point data arrays
    int type = scalarArray->GetDataType();
    vtkSmartPointer<vtkDataArray> values = vtkDataArray::CreateDataArray(type);
    values->SetNumberOfComponents(scalarArray->GetNumberOfComponents());
    values->SetName(scalarArray->GetName());
    outPD->AddArray(values);

    // additional data on ridge surface
    vtkSmartPointer<vtkDataArray> addData = 
        vtkDataArray::CreateDataArray(dataArray->GetDataType());
    addData->SetNumberOfComponents(dataArray->GetNumberOfComponents());
    addData->SetName(dataArray->GetName());
    outPD->AddArray(addData);


    // additional tet grid showing processed tetrahedra
    outgrid->ShallowCopy(input);
    vtkSmartPointer<vtkCellArray> tets = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkFloatArray> gradDotEVArray = vtkSmartPointer<vtkFloatArray>::New();
    gradDotEVArray->SetNumberOfComponents(1);
    gradDotEVArray->SetNumberOfTuples(numPts);
    gradDotEVArray->SetName("GradDotEV");
    outgrid->GetPointData()->AddArray(gradDotEVArray);

    vtkSmartPointer<vtkIntArray> triIndexArray = vtkIntArray::New();
    triIndexArray->SetNumberOfComponents(1);
    triIndexArray->SetName("TriangleIndex");
    outgrid->GetCellData()->AddArray(triIndexArray);

    QElapsedTimer timer;
    timer.start();

    // iterate over all cells, insert triangles forming the iso surface
    for( cellId=0; cellId<numCells; cellId++) 
    {
        vtkCell* cell = input->GetCell(cellId);
        if( cell && cell->GetCellType() == VTK_TETRA ) 
        {
            vtkSmartPointer<vtkIdList> pointIds = cell->GetPointIds();

            int ptIds[4];
            ptIds[0] = pointIds->GetId(0);
            ptIds[1] = pointIds->GetId(1);
            ptIds[2] = pointIds->GetId(2);
            ptIds[3] = pointIds->GetId(3);

            double ev[3*4], lev[3];
            eV->GetTuple(ptIds[0], &ev[0*3]);
            eV->GetTuple(ptIds[1], &ev[1*3]);
            eV->GetTuple(ptIds[2], &ev[2*3]);
            eV->GetTuple(ptIds[3], &ev[3*3]);
            getLargestEV(ev, lev);

            double gradDotEV[4];
            for (int n=0; n<4; ++n) 
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

            gradDotEVArray->SetValue(ptIds[0], gradDotEV[0]);
            gradDotEVArray->SetValue(ptIds[1], gradDotEV[1]);
            gradDotEVArray->SetValue(ptIds[2], gradDotEV[2]);
            gradDotEVArray->SetValue(ptIds[3], gradDotEV[3]);

            int triindex = 0;

            XYZ p0, p1, p2, p3;
            input->GetPoint(ptIds[0], &p0.x);
            input->GetPoint(ptIds[1], &p1.x);
            input->GetPoint(ptIds[2], &p2.x);
            input->GetPoint(ptIds[3], &p3.x);

            double v[4];
            scalarArray->GetTuple(ptIds[0], &v[0]);
            scalarArray->GetTuple(ptIds[1], &v[1]);
            scalarArray->GetTuple(ptIds[2], &v[2]);
            scalarArray->GetTuple(ptIds[3], &v[3]);

            double ew[4];
            eW->GetTuple(ptIds[0], &ew[0]);
            eW->GetTuple(ptIds[1], &ew[1]);
            eW->GetTuple(ptIds[2], &ew[2]);
            eW->GetTuple(ptIds[3], &ew[3]);

            double data[4][9];
            dataArray->GetTuple(ptIds[0], &data[0][0]);
            dataArray->GetTuple(ptIds[1], &data[1][0]);
            dataArray->GetTuple(ptIds[2], &data[2][0]);
            dataArray->GetTuple(ptIds[3], &data[3][0]);

            if( gradDotEV[0] < 0 ) triindex |= 1;
            if( gradDotEV[1] < 0 ) triindex |= 2;
            if( gradDotEV[2] < 0 ) triindex |= 4;
            if( gradDotEV[3] < 0 ) triindex |= 8;

            XYZ np0, np1, np2;
            double val0, val1, val2;
            double ew0, ew1, ew2;
            double d0[9], d1[9], d2[9];
            double t = this->EigenValueThreshold;
            bool triangle = false;

            switch( triindex )
            {
            case 0x00:
            case 0x0F:
                break;
            case 0x0E:
            case 0x01:                
                np0 = VertexInterpolate(0.0,p0,p1,gradDotEV[0],gradDotEV[1], v[0], v[1], val0, ew[0], ew[1], ew0, data[0], data[1], d0);
                np1 = VertexInterpolate(0.0,p0,p2,gradDotEV[0],gradDotEV[2], v[0], v[2], val1, ew[0], ew[2], ew1, data[0], data[2], d1);
                np2 = VertexInterpolate(0.0,p0,p3,gradDotEV[0],gradDotEV[3], v[0], v[3], val2, ew[0], ew[3], ew2, data[0], data[3], d2);
                if( applyFilter(val0, val1, val2, ew0, t, ew1, ew2) ) 
                {
                    InsertTriangle(np0, np1, np2, newPoints, newTris);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                    values->InsertNextTuple1(val2);
                    addData->InsertNextTuple(d0);
                    addData->InsertNextTuple(d1);
                    addData->InsertNextTuple(d2);
                    triangle = true;
                }
                break;
            case 0x0D:
            case 0x02:
                np0 = VertexInterpolate(0.0,p1,p0,gradDotEV[1],gradDotEV[0], v[1], v[0], val0, ew[1], ew[0], ew0, data[1], data[0], d0);
                np1 = VertexInterpolate(0.0,p1,p3,gradDotEV[1],gradDotEV[3], v[1], v[3], val1, ew[1], ew[3], ew1, data[1], data[3], d1);
                np2 = VertexInterpolate(0.0,p1,p2,gradDotEV[1],gradDotEV[2], v[1], v[2], val2, ew[1], ew[2], ew2, data[1], data[2], d2);
                if( applyFilter(val0, val1, val2, ew0, t, ew1, ew2) ) 
                {
                    InsertTriangle(np0, np1, np2, newPoints, newTris);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                    values->InsertNextTuple1(val2);
                    addData->InsertNextTuple(d0);
                    addData->InsertNextTuple(d1);
                    addData->InsertNextTuple(d2);
                    triangle = true;
                }
                break;
            case 0x0C:
            case 0x03:
                np0 = VertexInterpolate(0.0,p0,p3,gradDotEV[0],gradDotEV[3], v[0], v[3], val0, ew[0], ew[3], ew0, data[0], data[3], d0);
                np1 = VertexInterpolate(0.0,p0,p2,gradDotEV[0],gradDotEV[2], v[0], v[2], val1, ew[0], ew[2], ew1, data[0], data[2], d1);
                np2 = VertexInterpolate(0.0,p1,p3,gradDotEV[1],gradDotEV[3], v[1], v[3], val2, ew[1], ew[3], ew2, data[1], data[3], d2);
                if( applyFilter(val0, val1, val2, ew0, t, ew1, ew2) ) 
                {
                    InsertTriangle(np0, np1, np2, newPoints, newTris);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                    values->InsertNextTuple1(val2);
                    addData->InsertNextTuple(d0);
                    addData->InsertNextTuple(d1);
                    addData->InsertNextTuple(d2);
                    triangle = true;
                }
                np0 = np2;
                np2 = np1;
                np1 = VertexInterpolate(0.0,p1,p2,gradDotEV[1],gradDotEV[2], v[1], v[2], val0, ew[1], ew[2], ew0, data[1], data[2], d0);
                if( applyFilter(val0, val1, val2, ew0, t, ew1, ew2) ) 
                {
                    InsertTriangle(np0, np1, np2, newPoints, newTris);
                    values->InsertNextTuple1(val2);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                    addData->InsertNextTuple(d2);
                    addData->InsertNextTuple(d0);
                    addData->InsertNextTuple(d1);
                    triangle = true;
                }
                break;
            case 0x0B:
            case 0x04:
                np0 = VertexInterpolate(0.0,p2,p0,gradDotEV[2],gradDotEV[0], v[2], v[0], val0, ew[2], ew[0], ew0, data[2], data[0], d0);
                np1 = VertexInterpolate(0.0,p2,p1,gradDotEV[2],gradDotEV[1], v[2], v[1], val1, ew[2], ew[1], ew1, data[2], data[1], d1);
                np2 = VertexInterpolate(0.0,p2,p3,gradDotEV[2],gradDotEV[3], v[2], v[3], val2, ew[2], ew[3], ew2, data[2], data[3], d2);
                if( applyFilter(val0, val1, val2, ew0, t, ew1, ew2) ) 
                {
                    InsertTriangle(np0, np1, np2, newPoints, newTris);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                    values->InsertNextTuple1(val2);
                    addData->InsertNextTuple(d0);
                    addData->InsertNextTuple(d1);
                    addData->InsertNextTuple(d2);
                    triangle = true;
                }
                break;
            case 0x0A:
            case 0x05:
                np0 = VertexInterpolate(0.0,p0,p1,gradDotEV[0],gradDotEV[1], v[0], v[1], val0, ew[0], ew[1], ew0, data[0], data[1], d0);
                np1 = VertexInterpolate(0.0,p2,p3,gradDotEV[2],gradDotEV[3], v[2], v[3], val1, ew[2], ew[3], ew1, data[2], data[3], d1);
                np2 = VertexInterpolate(0.0,p0,p3,gradDotEV[0],gradDotEV[3], v[0], v[3], val2, ew[0], ew[3], ew2, data[0], data[3], d2);
                if( applyFilter(val0, val1, val2, ew0, t, ew1, ew2) ) 
                {
                    InsertTriangle(np0, np1, np2, newPoints, newTris);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                    values->InsertNextTuple1(val2);                
                    addData->InsertNextTuple(d0);
                    addData->InsertNextTuple(d1);
                    addData->InsertNextTuple(d2);
                    triangle = true;
                }
                np0 = np0;
                np2 = np1;
                np1 = VertexInterpolate(0.0,p1,p2,gradDotEV[1],gradDotEV[2], v[1], v[2], val2, ew[1], ew[2], ew2, data[1], data[2], d2);
                if( applyFilter(val0, val1, val2, ew0, t, ew1, ew2) ) 
                {
                    InsertTriangle(np0, np1, np2, newPoints, newTris);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val2);                
                    values->InsertNextTuple1(val1);
                    addData->InsertNextTuple(d0);
                    addData->InsertNextTuple(d2);
                    addData->InsertNextTuple(d1);
                    triangle = true;
                }
                break;
            case 0x09:
            case 0x06:
                np0 = VertexInterpolate(0.0,p0,p1,gradDotEV[0],gradDotEV[1], v[0], v[1], val0, ew[0], ew[1], ew0, data[0], data[1], d0);
                np1 = VertexInterpolate(0.0,p1,p3,gradDotEV[1],gradDotEV[3], v[1], v[3], val1, ew[1], ew[3], ew1, data[1], data[3], d1);
                np2 = VertexInterpolate(0.0,p2,p3,gradDotEV[2],gradDotEV[3], v[2], v[3], val2, ew[2], ew[3], ew2, data[2], data[3], d2);
                if( applyFilter(val0, val1, val2, ew0, t, ew1, ew2) ) 
                {
                    InsertTriangle(np0, np1, np2, newPoints, newTris);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                    values->InsertNextTuple1(val2);
                    addData->InsertNextTuple(d0);
                    addData->InsertNextTuple(d1);
                    addData->InsertNextTuple(d2);
                    triangle = true;
                }
                np0 = np0;
                np1 = VertexInterpolate(0.0,p0,p2,gradDotEV[0],gradDotEV[2], v[0], v[2], val1, ew[0], ew[2], ew1, data[0], data[2], d1);
                np2 = np2;
                if( applyFilter(val0, val1, val2, ew0, t, ew1, ew2) ) 
                {
                    InsertTriangle(np0, np1, np2, newPoints, newTris);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                    values->InsertNextTuple1(val2);
                    addData->InsertNextTuple(d0);
                    addData->InsertNextTuple(d1);
                    addData->InsertNextTuple(d2);
                    triangle = true;
                }
                break;
            case 0x07:
            case 0x08:
                np0 = VertexInterpolate(0.0,p3,p0,gradDotEV[3],gradDotEV[0], v[3], v[0], val0, ew[3], ew[0], ew0, data[3], data[0], d0);
                np1 = VertexInterpolate(0.0,p3,p2,gradDotEV[3],gradDotEV[2], v[3], v[2], val1, ew[3], ew[2], ew1, data[3], data[2], d1);
                np2 = VertexInterpolate(0.0,p3,p1,gradDotEV[3],gradDotEV[1], v[3], v[1], val2, ew[3], ew[1], ew2, data[3], data[1], d2);
                if( applyFilter(val0, val1, val2, ew0, t, ew1, ew2) ) 
                {
                    InsertTriangle(np0, np1, np2, newPoints, newTris);
                    values->InsertNextTuple1(val0);
                    values->InsertNextTuple1(val1);
                    values->InsertNextTuple1(val2);
                    addData->InsertNextTuple(d0);
                    addData->InsertNextTuple(d1);
                    addData->InsertNextTuple(d2);
                    triangle = true;
                }
                break;
            }
            if( triangle) {
                tets->InsertNextCell(cell);
                triIndexArray->InsertNextValue(triindex);
            }
        }

        this->UpdateProgress((float)cellId/numCells);
    }

    write_timer("TetRidge", "calculate", timer.elapsed());

    if( this->RegionThreshold > 0) 
    {
        timer.restart();

        std::vector<unsigned int> indices_merged;
        std::vector<XYZ> vertices_merged;
        std::vector<double> scalars_merged;
        vtkSmartPointer<vtkDataArray> data_merged = 
            vtkDataArray::CreateDataArray(dataArray->GetDataType());        
        data_merged->SetNumberOfComponents(dataArray->GetNumberOfComponents());

        mergeTriangles(newPoints, values, addData, 
            indices_merged, vertices_merged, scalars_merged, data_merged);

        std::vector<int> triLabels;
        std::vector<int> compSizes;
        labelComponents(&indices_merged, Manifold, &triLabels, &compSizes);

        std::cout << "Region triangle count threshold: " << RegionThreshold << std::endl;
        std::set<int> drawIds;

        for (int nComp = 0; nComp < compSizes.size(); ++nComp) {
            if (compSizes[nComp] > RegionThreshold)
                drawIds.insert(nComp);
        }

        std::cout << "indices_merged.size(): " << indices_merged.size() << std::endl;
        std::cout << "vertices_merged.size(): " << vertices_merged.size() << std::endl;
        std::cout << "scalars_merged.size(): " << scalars_merged.size() << std::endl;

        newPoints = vtkSmartPointer<vtkPoints>::New();
        newTris   = vtkSmartPointer<vtkCellArray>::New();
        values->SetNumberOfTuples(0);
        addData->SetNumberOfTuples(0);

        for (int t=0; t<triLabels.size(); ++t) {
            vtkIdType id[1];
            int vIdx; 
            vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();
            bool drawTri = (drawIds.find(triLabels[t]) != drawIds.end());
            if (drawTri) {

                vIdx = indices_merged[t*3+0];
                id[0] = newPoints->InsertNextPoint(&vertices_merged[vIdx].x);
                values->InsertNextTuple1(scalars_merged[vIdx]);
                addData->InsertNextTuple(data_merged->GetTuple(vIdx));
                tri->GetPointIds()->SetId(0, id[0]);

                vIdx = indices_merged[t*3+1];
                id[0] = newPoints->InsertNextPoint(&vertices_merged[vIdx].x);
                values->InsertNextTuple1(scalars_merged[vIdx]);
                addData->InsertNextTuple(data_merged->GetTuple(vIdx));
                tri->GetPointIds()->SetId(1, id[0]);

                vIdx = indices_merged[t*3+2];
                id[0] = newPoints->InsertNextPoint(&vertices_merged[vIdx].x);
                values->InsertNextTuple1(scalars_merged[vIdx]);
                addData->InsertNextTuple(data_merged->GetTuple(vIdx));
                tri->GetPointIds()->SetId(2, id[0]);

                newTris->InsertNextCell(tri);
            }
        }

        write_timer("TetRidge", "FilterBySize", timer.elapsed());

        data_merged->Delete();
    }

    values->Delete();
    addData->Delete();
    
    output->SetPoints(newPoints);
    output->SetPolys(newTris);
    output->Squeeze();

    outgrid->SetCells(VTK_TETRA, tets);
    outgrid->Squeeze();

    return 1;
}



bool vtkTetRidge::applyFilter(double val0, double val1, double val2, double ew0, double t, double ew1, double ew2)
{
    if( this->StrictFilter )
        return val0 >= this->MinDataValue && val0 <= this->MaxDataValue &&
               val1 >= this->MinDataValue && val1 <= this->MaxDataValue &&
               val2 >= this->MinDataValue && val2 <= this->MaxDataValue &&
               ew0 <= t && ew1 <= t && ew2 <= t;
    else
        return (val0 >= this->MinDataValue && val0 <= this->MaxDataValue) ||
               (val1 >= this->MinDataValue && val1 <= this->MaxDataValue) ||
               (val2 >= this->MinDataValue && val2 <= this->MaxDataValue) ||
               ew0 <= t || ew1 <= t || ew2 <= t;
}

XYZ vtkTetRidge::VertexInterpolate( double isovalue, XYZ p1, XYZ p2, 
                                    double valp1, double valp2, 
                                    double data1, double data2, double& val, 
                                    double ew1, double ew2, double& ew,
                                    double d1[9], double d2[9], double d[9])
{
    bool swp = ( p1.x < p2.x ) ||
                ( p1.x == p2.x && p1.y < p2.y) ||
                ( p1.y == p2.y && p1.z < p2.z);
    
    if( swp ) {
        swap<XYZ>(p1, p2);
        swap<double>(valp1, valp2);
        swap<double>(data1, data2);
        swap<double>(ew1, ew2);

        for( int i=0; i<9; i++ ){
            double t = d1[i];
            d1[i] = d2[i];
            d2[i] = t;
        }            
    }

    XYZ p;
    if( fabs(isovalue-valp1) < 0.00001) {
        p = p1;
        val = data1;
        ew = ew1;
        for( int i=0; i<9; i++ ){
            d[i] = d1[i];
        }            
    } else if (fabs(isovalue-valp2) < 0.00001) {
        p = p2;
        val = data2;
        ew = ew2;
        for( int i=0; i<9; i++ ){
            d[i] = d2[i];
        }            
    } else if (fabs(valp1-valp2) < 0.00001) {
        p = p1;
        val = data1;
        ew = ew1;
        for( int i=0; i<9; i++ ){
            d[i] = d1[i];
        }            
    } else {
        double mu = (isovalue - valp1) / (valp2-valp1);
        p.x = p1.x  + mu * (p2.x  - p1.x);
        p.y = p1.y  + mu * (p2.y  - p1.y);
        p.z = p1.z  + mu * (p2.z  - p1.z);
        val = data1 + mu * (data2 - data1);
        ew  = ew1   + mu * (ew2   - ew1);
        for( int i=0; i<9; i++ ){
            d[i] = d1[i] + mu * (d2[i] - d1[i]);
        }            
    }

    return p;
}

void vtkTetRidge::InsertTriangle(XYZ np0, XYZ np1, XYZ np2, vtkPoints * newPoints, vtkCellArray * newTris)
{
    newTris->InsertNextCell(3);
    newTris->InsertCellPoint( newPoints->InsertNextPoint(&np0.x));
    newTris->InsertCellPoint( newPoints->InsertNextPoint(&np1.x));
    newTris->InsertCellPoint( newPoints->InsertNextPoint(&np2.x));    
}

int vtkTetRidge::FillInputPortInformation(int port, vtkInformation* info)
{
    if ( port == 0 ) {
        info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid" );
        return 1;
    }
    return 0;
}

int vtkTetRidge::FillOutputPortInformation(int port, vtkInformation* info)
{
    if( port == 0 ){
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    } else if ( port == 1 ){
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    } 
    return 1;
}

void vtkTetRidge::getLargestEV(double *ev, double lev[3])
{
    // calculate largest eigen vector of cell
    double evSet[3][4*2];
    for (int n=0; n<4; ++n) {
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
            for(int n=0; n<4; ++n){
                sum += evSet[i][n] * evSet[j][n];
            }
            C[j][i] = sum/16.0;
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

void vtkTetRidge::mergeTriangles( vtkPoints* points, 
    vtkDataArray* values, 
    vtkDataArray* data_in, 
    std::vector<unsigned int>& indices, 
    std::vector<XYZ>& vertices, 
    std::vector<double>& scalars,
    vtkDataArray* data_merged )
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
        FKeys3 key(pos.x, pos.y, pos.z/*, this->MergeTolerance*/);

        // find vertex
        if (vertexInfo.find(key) == vertexInfo.end()) {
            // not found

            vertexInfo[key] = vertexID;

            // insert vertex and its corresponding data
            vertices.push_back(pos);
            scalars.push_back(values->GetTuple1(ptId));
            data_merged->InsertNextTuple(data_in->GetTuple(ptId));

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

void vtkTetRidge::labelComponents(std::vector<unsigned int>* indices, 
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

int vtkTetRidge::labelComponents_rek(std::vector<unsigned int>* indices,
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