/*=========================================================================

Program:   Visualization Toolkit
Module:    PowerCrust.h

Copyright (c) Michael Bußler

=========================================================================*/

#include "PowerCrust.h"
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

#include "vtkPowerCrustSurfaceReconstruction.h"

vtkStandardNewMacro(PowerCrust);

//==============================================================================
PowerCrust::PowerCrust()
{
    this->SetNumberOfOutputPorts(2);
}

//==============================================================================
PowerCrust::~PowerCrust()
{
    // clean up
}

//==============================================================================
int PowerCrust::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *outInfo2 = outputVector->GetInformationObject(1);

    // get the input and output
    vtkDataSet *input = vtkDataSet::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output2 = vtkPolyData::SafeDownCast(
        outInfo2->Get(vtkDataObject::DATA_OBJECT()));

    vtkIdType estimatedSize, numCells=input->GetNumberOfCells();
    vtkIdType numPts=input->GetNumberOfPoints();
    //vtkPoints *inPts=input->GetPoints();
    int numberOfPoints;
    vtkPointData *inPD=input->GetPointData(), *outPD = output->GetPointData();
    vtkCellData *inCD=input->GetCellData(), *outCD = output->GetCellData();

    if ( numPts < 1 )
    {
        vtkDebugMacro(<<"No input data.");
        return 1;
    }

    vtkSmartPointer<vtkPowerCrustSurfaceReconstruction> SurfaceReconstructor = 
        vtkSmartPointer<vtkPowerCrustSurfaceReconstruction>::New();
    SurfaceReconstructor->SetInputData ( input );
    SurfaceReconstructor->Update();

    output->DeepCopy(SurfaceReconstructor->GetOutput());
    output2->DeepCopy(SurfaceReconstructor->GetMedialSurface());

    return 1;
}

//----------------------------------------------------------------------------
int PowerCrust::FillOutputPortInformation( int port, vtkInformation* info )
{
    if ( port == 0 ) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
        return 1;
    }
    if ( port == 1 ) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
        return 1;
    }

    return 0;
}

//==============================================================================
void PowerCrust::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}
