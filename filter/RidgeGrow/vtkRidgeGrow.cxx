
#include "vtkRidgeGrow.h"

#include "vtkCellData.h"
#include "vtkCompositeDataIterator.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkDistancePolyDataFilter.h"
#include "vtkHausdorffDistancePointSetFilter.h"
#include "vtkMath.h"
#include <set>
#include <map>
#include <algorithm>

/* performance measure */
#include "timer.h"
#include <QElapsedTimer>

vtkStandardNewMacro(vtkRidgeGrow);

vtkRidgeGrow::vtkRidgeGrow()
{
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(3);

    growMesh = vtkSmartPointer<vtkPolyData>::New();
    ClearCellGrowth = false;

    this->TargetDistanceMethod = vtkHausdorffDistancePointSetFilter::POINT_TO_CELL;
}

vtkRidgeGrow::~vtkRidgeGrow()
{
}

int vtkRidgeGrow::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    }
    if (port == 1) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    }
    return 1;
}

int vtkRidgeGrow::FillOutputPortInformation(int port, vtkInformation* info)
{
    if ( port == 0 ) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
        return 1;
    }
    if ( port == 1 ) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
        return 1;
    }
    if ( port == 2 ) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
        return 1;
    }
    return 0;
}

int vtkRidgeGrow::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector)
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *inInfoRef = inputVector[1]->GetInformationObject(0);
    vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
    vtkInformation *outInfo2 = outputVector->GetInformationObject(2);
    
    vtkPolyData *input = vtkPolyData::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *inputRef = vtkPolyData::SafeDownCast(
        inInfoRef->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo0->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *outputFragment = vtkPolyData::SafeDownCast(
        outInfo1->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *outputPrev = vtkPolyData::SafeDownCast(
        outInfo2->Get(vtkDataObject::DATA_OBJECT()));

    // get the requested update times
    double upTime = outInfo0->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    
    int     numTimes = 0;
    double *inTimes  = NULL;
    bool noTimesteps = false;
    if( inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {
        inTimes = inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        numTimes = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        noTimesteps = numTimes <= 1;
    } else {
        noTimesteps = true;
    }

    if( noTimesteps ) {
        vtkErrorMacro("No timesteps to grow.");
        return 0;
    }

    int numPtsRef = inputRef->GetNumberOfPoints();
    int numCellsRef = inputRef->GetNumberOfPolys();

    vtkPointData *inPD = inputRef->GetPointData();
    vtkPointData *outPD = growMesh->GetPointData();
    vtkCellData  *outCD = growMesh->GetCellData();
    int numArrays = inPD->GetNumberOfArrays();
    vtkSmartPointer<vtkDoubleArray> timeValues;
    vtkSmartPointer<vtkDoubleArray> growthRate;
    //vtkSmartPointer<vtkDoubleArray> offsetVectors;
    vtkSmartPointer<vtkIntArray> newCells;
    vtkSmartPointer<vtkIntArray> labels;
    vtkSmartPointer<vtkDoubleArray> cellGrowth;


    QElapsedTimer timer;
    timer.start();

    // init/reset current mesh and ref point vector
    int resetTimestep = std::max<int>(std::min<int>(this->ResetOnTimestep, numTimes-1), 0);
    if( refPtIds.size() == 0 || upTime <= inTimes[resetTimestep])
    {
        refPtIds.clear();
        refPtIds.resize(numPtsRef, -1);
        refCellIds.clear();
        refCellIds.resize(numCellsRef, false);
        isolatedPt.clear();
        isolatedPt.resize(numPtsRef, true);
        
        growMesh = vtkSmartPointer<vtkPolyData>::New();
        growMesh->SetPoints(vtkPoints::New());
        growMesh->SetPolys(vtkCellArray::New());
        outPD = growMesh->GetPointData();
        outCD = growMesh->GetCellData();

        // copy point data arrays
        for( int i=0; i<numArrays; i++ ){
            int type = inPD->GetArray(i)->GetDataType();
            vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
            arr->SetNumberOfComponents(inPD->GetArray(i)->GetNumberOfComponents());
            arr->SetName(inPD->GetArray(i)->GetName());
            outPD->AddArray(arr);
            arr->Delete();
        }
        
        timeValues = vtkSmartPointer<vtkDoubleArray>::New();
        timeValues->SetNumberOfComponents(1);
        timeValues->SetName("AddTime");
        outPD->AddArray(timeValues);

        growthRate = vtkSmartPointer<vtkDoubleArray>::New();
        growthRate->SetNumberOfComponents(1);
        growthRate->SetName("GrowthRate");
        outPD->AddArray(growthRate);

        //offsetVectors = vtkSmartPointer<vtkDoubleArray>::New();
        //offsetVectors->SetNumberOfComponents(3);
        //offsetVectors->SetName("OffsetVector");
        //outPD->AddArray(offsetVectors);

        // flag new triangles
        newCells = vtkSmartPointer<vtkIntArray>::New();
        newCells->SetName("NewCell");
        newCells->SetNumberOfComponents(1);
        outCD->AddArray(newCells);

        labels = vtkSmartPointer<vtkIntArray>::New();
        labels->SetName("labels");
        labels->SetNumberOfComponents(1);
        outCD->AddArray(labels);

        cellGrowth = vtkSmartPointer<vtkDoubleArray>::New();
        cellGrowth->SetName("CellGrowth");
        cellGrowth->SetNumberOfComponents(1);
        outCD->AddArray(cellGrowth);

    }
    timeValues = vtkDoubleArray::SafeDownCast(outPD->GetArray("AddTime"));
    //offsetVectors = vtkDoubleArray::SafeDownCast(outPD->GetArray("OffsetVector"));
    growthRate = vtkDoubleArray::SafeDownCast(outPD->GetArray("GrowthRate"));
    newCells = vtkIntArray::SafeDownCast(outCD->GetArray("NewCell"));

    // calculate per point distance of input mesh to reference mesh
    vtkSmartPointer<vtkHausdorffDistancePointSetFilter> distanceFilter =
        vtkSmartPointer<vtkHausdorffDistancePointSetFilter>::New();
    distanceFilter->SetTargetDistanceMethod( this->TargetDistanceMethod );

    distanceFilter->SetInputData(0, inputRef );
    distanceFilter->SetInputData(1, input );
    distanceFilter->Update();

    vtkPolyData* res = distanceFilter->GetOutput(0);
    vtkDataArray* distArray = res->GetPointData()->GetArray("Distance");
    vtkIdTypeArray* nearestIdArray = vtkIdTypeArray::SafeDownCast(
        res->GetPointData()->GetArray("NearestPointId"));

    // copy active part of reference mesh to output
    vtkSmartPointer<vtkPoints> points = growMesh->GetPoints();
    vtkSmartPointer<vtkPoints> pointsFrag = vtkSmartPointer<vtkPoints>::New();
    
    std::vector<int> refPtIdsFrag;
    refPtIdsFrag.resize(numPtsRef, -1);

    // store current state of grow mesh
    vtkSmartPointer<vtkPoints> prevPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> prevCells = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> prevMesh = vtkSmartPointer<vtkPolyData>::New();
    prevPoints->DeepCopy(growMesh->GetPoints());
    prevCells->DeepCopy(growMesh->GetPolys());
    prevMesh->SetPoints(prevPoints);
    prevMesh->SetPolys(prevCells);
    prevMesh->BuildCells();

    outputPrev->DeepCopy(growMesh);

    vtkIdType firstNewPointId = -1;

    // compare input and ref point-wise and activate nearest points from ref
    if( distArray )
    {
        for( vtkIdType ptId=0; ptId<res->GetNumberOfPoints(); ptId++)
        {
            double dist = distArray->GetTuple1(ptId);
            
            // check if point is near point in ref mesh
            if( dist <= this->MinDistance ) 
            {
                double pos[3];
                inputRef->GetPoint(ptId, pos);

                //double pos2[3];
                //input->GetPoint(nearestIdArray->GetValue(ptId), pos2);

                //double offsetVector[3];
                //vtkMath::Subtract(pos2, pos, offsetVector);

                if( refPtIds[ptId] < 0 ) // point not inserted yet
                {
                    // insert point into output
                    vtkIdType newId = points->InsertNextPoint(pos);
                    refPtIds[ptId] = newId;

                    vtkIdType newIdFrag = pointsFrag->InsertNextPoint(pos);
                    refPtIdsFrag[ptId] = newIdFrag;

                    // copy point data
                    for( int i=0; i<numArrays; i++ ){
                        double data[9];
                        inPD->GetArray(i)->GetTuple(ptId, data);
                        outPD->GetArray(i)->InsertNextTuple(data);
                    }
                    timeValues->InsertNextValue(upTime);
                    //offsetVectors->InsertNextTuple(offsetVector);

                    if( firstNewPointId < 0)
                        firstNewPointId = newId; // new points are always appended, this is the first one
                }
                //else 
                //{
                //    offsetVectors->SetTuple(refPtIds[ptId], offsetVector);
                //}
            }
        }
    }

    outputFragment->SetPoints(pointsFrag);

    growthRate->SetNumberOfTuples(growMesh->GetNumberOfPoints());
    for( vtkIdType ptId=0; ptId<growMesh->GetNumberOfPoints(); ptId++)
        growthRate->SetValue(ptId, 0.0);

    if( prevMesh->GetNumberOfCells() > 0) {

        // estimate growth rate by calculating distance between prev and new grow mesh
        distanceFilter = vtkSmartPointer<vtkHausdorffDistancePointSetFilter>::New();
        distanceFilter->SetTargetDistanceMethod( this->TargetDistanceMethod );
        distanceFilter->SetInputData(0, growMesh );
        distanceFilter->SetInputData(1, prevMesh );
        distanceFilter->Update();

        res = distanceFilter->GetOutput(0);
        distArray = res->GetPointData()->GetArray("Distance");
   
        // copy point-wise distance to newly generated mesh
        if( distArray && firstNewPointId >= 0 )
        {
            //growthRate->SetNumberOfTuples(0);
            for( vtkIdType ptId=firstNewPointId; ptId<growMesh->GetNumberOfPoints(); ptId++)
            {
                double dist = distArray->GetTuple1(ptId);
                //if( dist <= this->MaxGrowthRate)
                    growthRate->SetValue(ptId, dist);
            }
        }
    }

    // Hello Paraview, something has changed here!!!
    timeValues->Modified();
    //offsetVectors->Modified();
    growthRate->Modified();

    // copy cells
    vtkSmartPointer<vtkCellArray> cells = growMesh->GetPolys();
    vtkCellArray* tris_ref = inputRef->GetPolys();
    tris_ref->InitTraversal();
    vtkSmartPointer<vtkIdList> cellPts = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> newCellPtIds, newCellPtIdsFrag;

    // reset previously inserted cells
    for(vtkIdType i=0; i < cells->GetNumberOfCells(); i++)
        newCells->SetValue(i, 0);
    
    // collect new cells
    std::vector<vtkIdType> newCellIds;

    vtkSmartPointer<vtkCellArray> cellsFrag = vtkSmartPointer<vtkCellArray>::New();

    cellGrowth = vtkDoubleArray::SafeDownCast(outCD->GetArray("CellGrowth"));

    vtkIdType cellId = 0;
    while( tris_ref->GetNextCell(cellPts) ) 
    {
        if( !refCellIds[cellId] ) // skip if cell was already inserted
        {
            // check if all cell points are active
            newCellPtIds = vtkSmartPointer<vtkIdList>::New();
            newCellPtIdsFrag = vtkSmartPointer<vtkIdList>::New();
            bool active = true;
            bool activeFrag = true;
            for( int i=0; i<cellPts->GetNumberOfIds(); i++){
                int newPtId = refPtIds[cellPts->GetId(i)];
                newCellPtIds->InsertNextId( newPtId );
                active &= (newPtId>=0); // active if all cell points are valid

                newPtId = refPtIdsFrag[cellPts->GetId(i)];
                newCellPtIdsFrag->InsertNextId( newPtId );
                activeFrag &= (newPtId>=0); // active if all cell points are valid
            }

            if( active ){
                vtkIdType newCellId = cells->InsertNextCell( newCellPtIds );
                newCells->InsertNextValue((int)1);
                refCellIds[cellId] = true;
                newCellIds.push_back(newCellId);
                cellGrowth->InsertNextValue(0);
            }
            if(activeFrag) {
                cellsFrag->InsertNextCell(newCellPtIdsFrag);
            }
        }
        cellId++;
    }

    outputFragment->SetPolys(cellsFrag);

    newCells->Modified();
    growMesh->BuildCells();
    growMesh->BuildLinks();

    // connected components of new cells
    int numCells = cells->GetNumberOfCells();
    int numNewCells = newCellIds.size();
    std::vector<bool> visited;
    visited.resize(numNewCells, false);

    labels = vtkIntArray::SafeDownCast(outCD->GetArray("labels"));
    labels->SetNumberOfValues(numCells);

    for( vtkIdType cellId=0; cellId<numCells; cellId++ )
    {
        labels->SetValue(cellId, 0);
        if( this->ClearCellGrowth)
            cellGrowth->SetValue(cellId, 0.0);
    }

    std::vector<std::vector<vtkIdType> > components;

    if( numCells != numNewCells )
    {
        int startCell = 0;
        int currentLabel = 1;
        std::vector<vtkIdType> stack;

        while( startCell < numNewCells)
        {
            stack.push_back(startCell);
            visited[startCell] = true;
            labels->SetValue(newCellIds[startCell], currentLabel);
            std::vector<vtkIdType> component;
            component.push_back(newCellIds[startCell]);

            while( !stack.empty())
            {
                int currentCell = stack.back();
                stack.pop_back();

                std::vector<vtkIdType> cellNeighbors;
                getCellNeighbors( growMesh, newCellIds[currentCell], cellNeighbors );

                std::vector<vtkIdType>::iterator cellNeighborsIter;
                for( cellNeighborsIter = cellNeighbors.begin(); cellNeighborsIter != cellNeighbors.end(); cellNeighborsIter++)
                {
                    std::vector<vtkIdType>::iterator cellIdIter;
                    cellIdIter = std::find(newCellIds.begin(), newCellIds.end(), *cellNeighborsIter );
                    if( cellIdIter != newCellIds.end())
                    {
                        int pos = cellIdIter - newCellIds.begin();
                        if( !visited[pos]) {
                            stack.push_back(pos);
                            labels->SetValue(*cellNeighborsIter, currentLabel);
                            component.push_back(*cellNeighborsIter);
                            visited[pos] = true;
                        }
                    }
                }
            }
            while( startCell < numNewCells && visited[startCell]) {
                startCell++;
            }
            currentLabel++;
            components.push_back(component);
        }
    }
    labels->Modified();

    vtkSmartPointer<vtkIntArray> connArray = vtkSmartPointer<vtkIntArray>::New();
    connArray->SetName("conn");
    connArray->SetNumberOfComponents(1);
    connArray->SetNumberOfTuples(growMesh->GetNumberOfPoints());
    for( int i=0; i<growMesh->GetNumberOfPoints(); i++)
        connArray->SetValue(i,0);

    // estimate maximum growth rate per component
    std::vector<vtkIdType> nonIsoPoints;
    for( std::vector<std::vector<vtkIdType> >::iterator componentsIter = components.begin();
         componentsIter != components.end(); componentsIter++)
    {
        double rate = 0.0;
        bool unconnected = true;
        for( std::vector<vtkIdType>::iterator componentIter = componentsIter->begin();
            componentIter != componentsIter->end(); componentIter++)
        {
            vtkSmartPointer<vtkIdList> cellPointIds =
                vtkSmartPointer<vtkIdList>::New();
            growMesh->GetCellPoints(*componentIter, cellPointIds);
            for(vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++)
            {
                vtkIdType ptId = cellPointIds->GetId(i);
                double pointRate = growthRate->GetValue(ptId);
                rate = std::max( rate, pointRate);
                if( ptId < firstNewPointId && !isolatedPt[ptId]){
                    unconnected = false;
                    connArray->SetValue(ptId, 1);
                }
                nonIsoPoints.push_back(ptId);
            }
        }
        // write back cell growth rate per component
        for( std::vector<vtkIdType>::iterator componentIter = componentsIter->begin();
            componentIter != componentsIter->end(); componentIter++)
        {
            if( unconnected )
                cellGrowth->SetValue(*componentIter, 0.0);
            else 
                cellGrowth->SetValue(*componentIter, rate);
        }

        //int pos = componentsIter - components.begin();
        //std::cout<< "component: " << pos <<" rate: "<<rate <<std::endl;
    }

    for( size_t i=0; i<nonIsoPoints.size(); i++)
        isolatedPt[nonIsoPoints[i]] = false;


    write_timer("ridgeGrow", "calculate", timer.elapsed());

    cellGrowth->Modified();

    output->ShallowCopy(growMesh);
    output->Squeeze();

    vtkSmartPointer<vtkIntArray> isolatedArray = vtkSmartPointer<vtkIntArray>::New();
    isolatedArray->SetName("isolated");
    isolatedArray->SetNumberOfComponents(1);
    isolatedArray->SetNumberOfTuples(output->GetNumberOfPoints());
    for( int i=0; i<output->GetNumberOfPoints(); i++)
        isolatedArray->SetValue(i,isolatedPt[i]);

    output->GetPointData()->AddArray(isolatedArray);
    output->GetPointData()->AddArray(connArray);

    return 1;
}

void vtkRidgeGrow::getCellNeighbors( vtkPolyData* mesh, vtkIdType cellId, std::vector<vtkIdType>& cellNeighbors)
{    
    vtkSmartPointer<vtkIdList> cellPointIds =
        vtkSmartPointer<vtkIdList>::New();
    mesh->GetCellPoints(cellId, cellPointIds);

    for(vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++)
    {
        vtkSmartPointer<vtkIdList> idList =
            vtkSmartPointer<vtkIdList>::New();
        idList->InsertNextId(cellPointIds->GetId(i));

        //add the other edge point
        if(i+1 == cellPointIds->GetNumberOfIds())
        {
            idList->InsertNextId(cellPointIds->GetId(0));
        }
        else
        {
            idList->InsertNextId(cellPointIds->GetId(i+1));
        }

        //get the neighbors of the cell
        vtkSmartPointer<vtkIdList> neighborCellIds =
            vtkSmartPointer<vtkIdList>::New();

        mesh->GetCellNeighbors(cellId, idList, neighborCellIds);

        for(vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); j++)
        {
            cellNeighbors.push_back( neighborCellIds->GetId(j) );
        }
    }

}
