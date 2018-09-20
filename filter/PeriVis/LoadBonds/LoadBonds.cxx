/*=========================================================================

Program:   Visualization Toolkit
Module:    LoadBonds.h

Copyright (c) Michael Bußler

=========================================================================*/

#include "LoadBonds.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"
#include "vtkSmartPointer.h"
#include "pqFileDialog.h"

#include <string>
#include <sstream>

vtkStandardNewMacro(LoadBonds);

//==============================================================================
LoadBonds::LoadBonds()
{
    this->CalculateBondLength = true;
    this->CalculateBrokenBonds = false;
    this->CalculateRefBondLength = true;
    this->FilterBondsByLength = true;
    this->MaxBondLength = 1.0;

    this->m_requested_time = 0;

    this->NumberOfTimeSteps = 0;
    this->m_initialMesh = NULL;

    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::SCALARS);
}

//==============================================================================
LoadBonds::~LoadBonds()
{
    // clean up
}

int LoadBonds::RequestInformation(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
    if( m_FileNames.empty())
    {
        // show file open dialog so that the user can select the bond file(s)
        pqFileDialog browseBondFileDialog( NULL, NULL, 
            "Select Bonds file(s)",
            ".", "Bonds file (*.csv)" );
        browseBondFileDialog.setFileMode(pqFileDialog::ExistingFile);
        if(pqFileDialog::Accepted == browseBondFileDialog.exec()) 
        { 
            QStringList files = browseBondFileDialog.getSelectedFiles();
            if(!files.empty()) {
                RemoveAllFileNames();
                for( int i=0; i<files.size(); i++)
                    AddFileName( files[i].toStdString().c_str() );
            }
        }
    }

    ParseBondsFile();

    //if( NumberOfTimeSteps > 1) 
    //{
    //    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    //    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), 
    //        &this->TimestepValues[0], this->NumberOfTimeSteps);

    //    double timeRange[2] = { 
    //        TimestepValues[0], 
    //        TimestepValues[NumberOfTimeSteps-1]
    //    };

    //    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
    //        timeRange, 2);

    //}

    return 1;
}

void LoadBonds::ParseBondsFile()
{
    TimestepValues.clear();
    m_timesteps.clear();

    for( int i=0; i<m_FileNames.size(); i++)
    {
        Timestep ts;
        ts.time = -1;
        ts.fileId = -1;
    
        ifstream inFile;
        inFile.open(m_FileNames[i].c_str(), ifstream::in);

        std::string line;
        std::getline(inFile, line);
        if( !line.empty() && line[0] == '#' ) 
        {
            std::vector<std::string> tokens;
            tokenize(line, tokens);
            if( tokens.size() == 6)
            {
                ts.time = stod(tokens[5]);
                ts.fileId = i;
            }
        }

        inFile.close();

        if( ts.time >= 0)
            m_timesteps.push_back(ts);
    }

    NumberOfTimeSteps = m_timesteps.size();

    for( int i=0; i<NumberOfTimeSteps; i++)
        TimestepValues.push_back(m_timesteps[i].time);
}

int LoadBonds::FindClosestTimeStep(double requestedTimeValue)
{
    int ts = 0;
    double mindist = std::abs(TimestepValues[0] - requestedTimeValue);

    for (int i = 0; i < NumberOfTimeSteps; i++) {

        double tsv = TimestepValues[i];
        double dist = std::abs(tsv - requestedTimeValue);
        if (dist < mindist) {
            mindist = dist;
            ts = i;
        }
    }
    return ts;
}

//==============================================================================
int LoadBonds::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkDataSet *input = vtkDataSet::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkPoints *newPoints;
    vtkCellArray *newLines;
    vtkIdType estimatedSize, numCells=input->GetNumberOfCells();
    vtkIdType numPts=input->GetNumberOfPoints();
    //vtkPoints *inPts=input->GetPoints();
    int numberOfPoints;
    vtkPointData *inPD=input->GetPointData(), *outPD = output->GetPointData();
    vtkCellData *inCD=input->GetCellData(), *outCD = output->GetCellData();

    vtkIdTypeArray *pointIds = vtkIdTypeArray::SafeDownCast(
        this->GetInputArrayToProcess(0, inputVector));

    if ( numPts < 1 || pointIds == NULL )
    {
        vtkDebugMacro(<<"No data to convert");
        return 1;
    }

    // copy point data arrays
    int numArrays = inPD->GetNumberOfArrays();
    for( int i=0; i<numArrays; i++ ){
        int type = inPD->GetArray(i)->GetDataType();
        vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
        arr->SetNumberOfComponents(inPD->GetArray(i)->GetNumberOfComponents());
        arr->SetName(inPD->GetArray(i)->GetName());
        outPD->AddArray(arr);
        arr->Delete();
    }

    // create data arrays for vertex cell data
    int numVertexArrays = inCD->GetNumberOfArrays();
    for( int i=0; i<numVertexArrays; i++ ){
        int type = inCD->GetArray(i)->GetDataType();
        vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
        arr->SetNumberOfComponents(inCD->GetArray(i)->GetNumberOfComponents());
        arr->SetName(inCD->GetArray(i)->GetName());
        outPD->AddArray(arr);
        arr->Delete();
    }
    
    // convert to polydata
    newPoints = vtkPoints::New();
        
    for( vtkIdType ptId=0; ptId<numPts; ptId++) 
    {
        // get ref point id from data array
        //vtkIdType refPtId = pointIds->LookupValue(ptId+1); // shift by one
        //newPoints->InsertNextPoint( input->GetPoint(refPtId));

        newPoints->InsertNextPoint(input->GetPoint(ptId));

        // copy point data
        for( int i=0; i<numArrays; i++ ){
            double data[9];
            inPD->GetArray(i)->GetTuple(ptId, data);
            outPD->GetArray(i)->InsertNextTuple(data);
        }

        // copy vertex data to points
        for( int i=0; i<numVertexArrays; i++ ){
            double data[9];
            inCD->GetArray(i)->GetTuple(ptId, data);
            outPD->GetArray(numArrays+i)->InsertNextTuple(data);
        }

        this->UpdateProgress((0.25*ptId)/numPts);
    }

    output->SetPoints(newPoints);

    int timestepToLoad = 0;
    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())
        && NumberOfTimeSteps > 0 ) {

            double requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

            timestepToLoad = FindClosestTimeStep(requestedTimeValue);
            output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), TimestepValues[timestepToLoad]);

            std::cout << "Loading timestep " << timestepToLoad << " with t=" << TimestepValues[timestepToLoad] << std::endl;
    }

    // currently only lines (bonds) are supported
    newLines = vtkCellArray::New();

    vtkIdTypeArray *pointIdsNew = vtkIdTypeArray::SafeDownCast(
        outPD->GetArray(pointIds->GetName()));

    // TODO: use FindClosestTimeStep here...
    Timestep& ts = m_timesteps[timestepToLoad];

    ifstream inFile;
    inFile.open(m_FileNames[ts.fileId], ifstream::in);

    // skip first line
    std::string line;
    getline(inFile, line);

    while (inFile.good()) 
    {
        getline(inFile, line);
        std::vector<std::string> tokens;
        tokenize(line, tokens);

        if( tokens.size() > 0 ) 
        {
            int ptId1 = std::stoi(tokens[0]);
            vtkIdType id1 = pointIdsNew->LookupValue(ptId1+1);

            if( id1 > 0 ) 
            {
                for( int i=1; i<tokens.size(); i++ ) 
                {
                    int ptId2 = std::stoi(tokens[i]);
                    vtkIdType id2 = pointIdsNew->LookupValue(ptId2+1);

                    if( id2 > 0) {
                        newLines->InsertNextCell(2);
                        newLines->InsertCellPoint( id1 );
                        newLines->InsertCellPoint( id2 );
                    }
                }
            }
        }
    }

    output->SetLines(newLines);
    this->UpdateProgress(0.5);

    //
    // Update output and release memory
    //

    if( timestepToLoad == 0 && !m_initialMesh )
    {
        m_initialMesh = vtkSmartPointer<vtkPolyData>::New();
        m_initialMesh->SetPoints(newPoints);
        m_initialMesh->SetLines(newLines);
    }

    newPoints->Delete();
    newLines->Delete();
    
    ProcessBondData(timestepToLoad, output);

    output->Squeeze();

    return 1;
}

void LoadBonds::ProcessBondData(int timestepToLoad, vtkPolyData* pd)
{
    // bond active flag: 1:active; 0:broken
    vtkSmartPointer<vtkFloatArray> bondActive = vtkSmartPointer<vtkFloatArray>::New();;
    bondActive->SetName("Bond Active");
    bondActive->SetNumberOfComponents(1);

    // length of bond in current timestep
    vtkSmartPointer<vtkFloatArray> bondLength = vtkSmartPointer<vtkFloatArray>::New();
    bondLength->SetName("Bond Length");
    bondLength->SetNumberOfComponents(1);

    // length of bond in reference timestep
    vtkSmartPointer<vtkFloatArray> refBondLength = vtkSmartPointer<vtkFloatArray>::New();
    refBondLength->SetName("Reference Bond Length");
    refBondLength->SetNumberOfComponents(1);

    /* Re-insert broken bonds by comparing current to reference state */
    if( this->GetCalculateBrokenBonds() && m_initialMesh)
    {
        int numBondsRef = m_initialMesh->GetNumberOfLines();
        bondActive->SetNumberOfTuples(numBondsRef);

        if( timestepToLoad == 0)
        {
            for( int cellId=0; cellId<numBondsRef; cellId++ )
                bondActive->SetValue(cellId, 1.0); // at T=0 all bonds are active (1)
        }
        else
        {
            // initalize result array
            for( int cellId=0; cellId<numBondsRef; cellId++ ) {
                bondActive->SetValue(cellId, 0.0); // initially all bonds are non-active (1)
            }

            // use reference lines set with changes to current state
            vtkSmartPointer<vtkCellArray> reflines = vtkSmartPointer<vtkCellArray>::New();
            reflines->DeepCopy(this->m_initialMesh->GetLines()); // use bonds of initial timestep
            pd->BuildLinks(); // current state

            // bonds per point id in current and ref data set
            vtkSmartPointer<vtkIdList> cellIdsRef = vtkSmartPointer<vtkIdList>::New();
            vtkSmartPointer<vtkIdList> cellIdsCurr = vtkSmartPointer<vtkIdList>::New();

            // now we need to compare the vertices of each bond in current and ref config against each other
            // assume that there are no more particles and bonds in the current config than in ref

            int numPoints = pd->GetNumberOfPoints();
            for( vtkIdType ptId=0; ptId<numPoints; ptId++) {

                vtkSmartPointer<vtkIdList> bondsIds = vtkSmartPointer<vtkIdList>::New();
                vtkSmartPointer<vtkIdList> neighborsCurr = getConnectedPoints(pd, ptId);
                vtkSmartPointer<vtkIdList> neighborsRef = getConnectedPoints(m_initialMesh, ptId, bondsIds);

                // compare ids to ref config
                for(vtkIdType i = 0; i<neighborsRef->GetNumberOfIds(); i++)
                {
                    vtkIdType neighIdRef = neighborsRef->GetId(i);

                    // check if id is also in current neighbor set
                    for( vtkIdType j = 0; j<neighborsCurr->GetNumberOfIds(); j++)
                    {
                        vtkIdType neighIdCurr = neighborsCurr->GetId(j);
                        if( neighIdCurr == neighIdRef ) { // bond is still active
                            vtkIdType bondId = bondsIds->GetId(i);
                            bondActive->SetValue(bondId, 1.0);
                            break;
                        }
                    }
                }
                this->UpdateProgress (0.5+0.25*ptId/(float)numPoints);
            }
            pd->SetLines(reflines);
        }
        pd->GetCellData()->AddArray(bondActive);
    }

    /* Estimate length of bonds */
    if(this->CalculateBondLength) 
    {
        pd->BuildCells();
        size_t numBonds = pd->GetNumberOfCells();

        bondLength->SetNumberOfTuples(numBonds);

        for(vtkIdType bondId = 0; bondId < numBonds; bondId++)
        {
            vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
            pd->GetCellPoints( bondId, pointIdList);
            double p0[3], p1[3];
            /* now use point coordinates from current state */
            pd->GetPoint( pointIdList->GetId(0), p0);
            pd->GetPoint( pointIdList->GetId(1), p1);
            double d2 = vtkMath::Distance2BetweenPoints(p0,p1);
            bondLength->SetValue(bondId, sqrt(d2));
        
            this->UpdateProgress (0.75+0.1*bondId/(float)numBonds);
        }

        pd->GetCellData()->AddArray(bondLength);
    }

    /* Estimate length of bonds in reference state*/
    if( this->CalculateRefBondLength && m_initialMesh )
    {
        pd->BuildCells();
        size_t numBonds = pd->GetNumberOfCells();

        refBondLength->SetNumberOfTuples(numBonds);

        for( int bondId=0; bondId<numBonds; bondId++ ) 
        {
            vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
            pd->GetCellPoints( bondId, pointIdList);
            double p0[3], p1[3];
            /* now use point coordinates from reference state */
            m_initialMesh->GetPoint( pointIdList->GetId(0), p0);
            m_initialMesh->GetPoint( pointIdList->GetId(1), p1);
            double d2 = vtkMath::Distance2BetweenPoints(p0,p1);
            refBondLength->SetValue(bondId, sqrt(d2));
        
            this->UpdateProgress (0.85+0.1*bondId/(float)numBonds);
        }

        pd->GetCellData()->AddArray(refBondLength);
    }

    /* Filter bonds by certain bond length */
    if( this->FilterBondsByLength && this->CalculateRefBondLength && m_initialMesh )
    {
        // new bonds, filtered by length
        vtkSmartPointer<vtkCellArray> newLines = vtkSmartPointer<vtkCellArray>::New();

        // new cell data arrays
        vtkSmartPointer<vtkFloatArray> bondActive_new = vtkSmartPointer<vtkFloatArray>::New();;
        bondActive_new->SetName("Bond Active");
        bondActive_new->SetNumberOfComponents(1);
        vtkSmartPointer<vtkFloatArray> bondLength_new = vtkSmartPointer<vtkFloatArray>::New();
        bondLength_new->SetName("Bond Length");
        bondLength_new->SetNumberOfComponents(1);
        vtkSmartPointer<vtkFloatArray> refBondLength_new = vtkSmartPointer<vtkFloatArray>::New();
        refBondLength_new->SetName("Reference Bond Length");
        refBondLength_new->SetNumberOfComponents(1);

        // current lines
        vtkCellArray* lines = pd->GetLines();
        size_t numBonds = lines->GetNumberOfCells();
        vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();

        // actual filtering
        for(vtkIdType bondId = 0; bondId < numBonds; bondId++)
        {
            float length = refBondLength->GetValue(bondId);
            if( length <= this->MaxBondLength )
            {
                pd->GetCellPoints(bondId, pts);
                //lines->GetCell(bondId, pts);
                newLines->InsertNextCell(pts);

                // copy cell data
                if( CalculateBrokenBonds)
                    bondActive_new->InsertNextValue( bondActive->GetValue(bondId));
                if( CalculateBondLength )
                    bondLength_new->InsertNextValue( bondLength->GetValue(bondId));
                refBondLength_new->InsertNextValue( refBondLength->GetValue(bondId));
            }
        }

        // replace data with new data
        pd->SetLines(newLines);
        if( CalculateBrokenBonds ){
            pd->GetCellData()->RemoveArray("Bond Active");
            pd->GetCellData()->AddArray(bondActive_new);
        }
        if( CalculateBondLength ){
            pd->GetCellData()->RemoveArray("Bond Length");
            pd->GetCellData()->AddArray(bondLength_new);
        }
        pd->GetCellData()->RemoveArray("Reference Bond Length");
        pd->GetCellData()->AddArray(refBondLength_new);
    }


}

void LoadBonds::AddFileName(const char* fname)
{
    m_FileNames.push_back(fname);
}

void LoadBonds::RemoveAllFileNames()
{
    m_FileNames.clear();
}

void LoadBonds::tokenize(const std::string line, std::vector<std::string> &tokens)
{
    std::stringstream iss(line);
    std::string token;
    while( iss >> token)
        tokens.push_back(token);
}


vtkSmartPointer<vtkIdList> LoadBonds::getConnectedPoints( vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds /*= NULL*/)
{
    vtkSmartPointer<vtkIdList> cellIds  = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> neighPts = vtkSmartPointer<vtkIdList>::New();

    pd->GetPointCells( ptId, cellIds);

    for(vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
    {
        vtkIdType cellId = cellIds->GetId(i);

        vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
        pd->GetCellPoints( cellId, pointIdList);

        vtkIdType ptIdN = pointIdList->GetId(0);
        if( ptIdN == ptId)
            ptIdN = pointIdList->GetId(1);

        neighPts->InsertNextId(ptIdN);

        if( bondIds )
            bondIds->InsertNextId(cellId);
    }

    return neighPts;
}

int LoadBonds::RequestUpdateExtent(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector )
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    
    //if( !m_initialMesh )
    //    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->TimestepValues[0]);

    return 1;
}

//----------------------------------------------------------------------------
int LoadBonds::FillInputPortInformation(int, vtkInformation *info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
}

//==============================================================================
void LoadBonds::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}
