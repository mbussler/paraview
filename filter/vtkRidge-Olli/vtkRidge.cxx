
#include "vtkRidge.h"

#include "vtkObjectFactory.h" //for new() macro
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkImageData.h"
#include "isosurface.h"

#include "vtkFloatArray.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Creates isosurface                                                                       ///
/// Creates an isosurface with a threshold passed in ParaView 'Create isosurface at' variable and   ///
/// scalar data on GPU using Marching Cubes algorithm.                                              ///
/// @param[in] c Isosurface threshold value                                                         ///
/// @param[in] data Scalar data set                                                                 ///
/// @param[in] dim Dimensions of data set in x,y,z                                                  ///
/// @param[out] numVerts Adress of uint to store number of generated vertices                       ///
/// @param[out] vertices Adress of pointer to receive array of doubles containing vertices x,y,z    ///
///////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" bool cudaIsosurface(const double c, const double* data, double* grads, uint dim[3], uint* numVerts, double** vertices, double** gradients);

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Creates isosurface                                                                       ///
/// Creates an array of gradient vectors based on the scalar field using GPGPU                      ///
/// @param[in] data Scalar data set                                                                 ///
/// @param[in] dim Dimensions of data set in x,y,z                                                  ///
/// @param[out] gradients Adress of pointer to receive array of doubles containing gradients x,y,z  ///
///////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" bool cudaGradient(const double* data, uint dim[3], double** gradients);

vtkStandardNewMacro(vtkRidge);

//-----------------------------------------------------------------------------
vtkRidge::vtkRidge()
{
}

//-----------------------------------------------------------------------------
vtkRidge::~vtkRidge()
{
}

//----------------------------------------------------------------------------
int vtkRidge::FillInputPortInformation( int port, vtkInformation* info )
{
    if ( port == 0 ) {
        info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData" );
        return 1;
    }
    // if (!this->Superclass::FillInputPortInformation(port, info)) {
    // 	return 0;
    // }
    return 0;
}

//----------------------------------------------------------------------------
int vtkRidge::RequestData(vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    // get the input and output objects from info
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    //this following object needs to be extracted from the input
    vtkSmartPointer<vtkImageData> in_grid;

    //get input data
    if (!inInfo->Has(vtkDataObject::DATA_OBJECT())) {
        cout << "Input has no data object. No calculation done." << endl;
        return 1;
    }
    else {
        in_grid = vtkImageData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
        double* origin;
        origin = inInfo->Get(vtkDataObject::ORIGIN());
        int* extent;
        extent = inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
        this->dim[0] = extent[1]+1;
        this->dim[1] = extent[3]+1;
        this->dim[2] = extent[5]+1;
    }

    //------------------------------------------------------------------------------------------------------------------------
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkDataSet *output = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkCellArray> out_topology = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> out_points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray> out_data = vtkSmartPointer<vtkFloatArray>::New();

    vtkDataArray *array = this->GetInputArrayToProcess(0, inputVector);
    if( !array )
    {
        vtkErrorMacro("No input array selected");
        return 0;
    } else if( !array->GetDataType() == VTK_DOUBLE || array->GetNumberOfComponents() != 1) {
        vtkErrorMacro("A scalar array of type double is required!");
        return 0;
    }

    in_grid->GetPointData()->SetActiveScalars( array->GetName());

    double* in_data = static_cast<double *>(in_grid->GetScalarPointer());

    //initialize isosurface
    Isosurface test;
    memcpy(test.m_dim, this->dim, 3*sizeof(uint));
    test.m_data = in_data;
    test.m_grad_data = NULL;

    //duh
    test.generateGradients();
#ifndef ISO_SHOW_GRADIENTS
    test.generateTriangles(isoThreshold);
#endif

    //cleanup
    test.m_data = NULL;

    vtkSmartPointer<vtkPolyData> opd = test.GetPolyDataPointer();

    //send to pipeline
    output->ShallowCopy(opd);
    return 1;
}

////////// External Operators /////////////

void vtkRidge::PrintSelf(ostream &os, vtkIndent indent)
{
}
