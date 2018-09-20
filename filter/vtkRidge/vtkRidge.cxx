
#include "vtkRidge.h"



#include "stdfunc.h"
#include "isosurface.h"
#include <iomanip>

#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h" //for new() macro
#include "vtkPixel.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTriangle.h"

//#include <sys/time.h>

/*
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
extern bool cudaIsosurface(const double c, const double* in_data, double* in_grads, int in_size[3], uint* out_numVerts, double** out_verts, double** out_grads);

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Creates isosurface                                                                       ///
/// Creates an array of gradient vectors based on the scalar field using GPGPU                      ///
/// @param[in] data Scalar data set                                                                 ///
/// @param[in] dim Dimensions of data set in x,y,z                                                  ///
/// @param[out] gradients Adress of pointer to receive array of doubles containing gradients x,y,z  ///
///////////////////////////////////////////////////////////////////////////////////////////////////////
extern bool cudaGradient(const double* data, uint dim[3], double** gradients);

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Creates isosurface                                                                       ///
/// Creates an array of gradient vectors based on the scalar field using GPGPU                      ///
/// @param[in] data Scalar data set                                                                 ///
/// @param[in] dim Dimensions of data set in x,y,z                                                  ///
/// @param[out] gradients Adress of pointer to receive array of doubles containing gradients x,y,z  ///
///////////////////////////////////////////////////////////////////////////////////////////////////////
extern bool cudaRidge(uint dim[3], const double* data, double** gradients, uint* numVerts, double** vertices);
*/
vtkStandardNewMacro(vtkRidge)

//-----------------------------------------------------------------------------
vtkRidge::vtkRidge()
{
    this->SetNumberOfOutputPorts(2);

    // by default process active point scalars
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::SCALARS);

    // by default process active point vectors
    this->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::VECTORS);

    // by default process active point vectors
    this->SetInputArrayToProcess(2,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::TENSORS);

    // by default process active point scalars
    this->SetInputArrayToProcess(3,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::SCALARS);

    this->Manifold = true;
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
    return 0;
}

//----------------------------------------------------------------------------
int vtkRidge::FillOutputPortInformation( int port, vtkInformation* info )
{
    if ( port == 0 ) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
        return 1;
    }
    if ( port == 1 ) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData" );
        return 1;
    }
    return 0;
}

//void vtkRidge::SetIsoThreshold(double iT) {
//    isoThreshold = iT;
//}

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
    //origin = inInfo->Get(vtkDataObject::ORIGIN());
	int* extent;
	int dom_dim[3];
	in_grid->GetDimensions(dom_dim);
	cout << "dom_dim: " << dom_dim[0] << " " << dom_dim[1] << " " << dom_dim[2] << endl;
	// extent = inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT());
	// cout << "extent: " << extent[0] << " " << extent[1] << " " << extent[2] << endl;
	this->dim[0] = dom_dim[0];
	this->dim[1] = dom_dim[1];
	this->dim[2] = dom_dim[2];
    }

    //------------------------------------------------------------------------------------------------------------------------
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkDataSet *output = vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
    vtkDataSet *output1 = vtkDataSet::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));
    
//    vtkSmartPointer<vtkCellArray> out_topology = vtkSmartPointer<vtkCellArray>::New();
//    vtkSmartPointer<vtkPoints> out_points = vtkSmartPointer<vtkPoints>::New();
//    vtkSmartPointer<vtkFloatArray> out_data = vtkSmartPointer<vtkFloatArray>::New();

    vtkDataArray *scalarArray = this->GetInputArrayToProcess(0, inputVector);
    if( !scalarArray )
    {
        vtkErrorMacro("No input array selected");
        return 0;
    }

    double *in_scalar = 0, 
           *in_grad   = 0, 
           *in_hesse  = 0;
    bool allocated_scalar_data = false,
         allocated_vector_data = false,
         allocated_tensor_data = false;

    in_grid->GetPointData()->SetActiveScalars( scalarArray->GetName());
    getData(in_grid, in_scalar, 1, allocated_scalar_data );

    vtkDataArray *vectorArray = this->GetInputArrayToProcess(1,inputVector);
    if( vectorArray && this->UseInput ){
        in_grid->GetPointData()->SetActiveScalars( vectorArray->GetName());
        getData(in_grid, in_grad, 3, allocated_vector_data );

        if( this->Valley ){ // flip sign of gradient
            const size_t nValues = dim[0]*dim[1]*dim[2]*3;
            for( int i=0; i<nValues; i++)
                in_grad[i]=-in_grad[i];
        }
    }

    vtkDataArray *tensorArray = this->GetInputArrayToProcess(2,inputVector);
    if( tensorArray && this->UseInput ){
        in_grid->GetPointData()->SetActiveScalars( tensorArray->GetName());
        getData(in_grid, in_hesse, 9, allocated_tensor_data );
    }

    int extent[6];
    in_grid->GetExtent(extent);
    double origin[3];
    in_grid->GetOrigin(origin);
    double spacing[3];
    in_grid->GetSpacing(spacing);

    //initialize isosurface
    Isosurface test;

    test.SetOrigin(origin);
    test.SetSpacing(spacing);
    test.SetExtent(extent);
    memcpy(test.m_dim, this->dim, 3*sizeof(uint));

    test.m_data = in_scalar;
    test.m_grad_data = in_grad;
    test.m_hesse_data = in_hesse;
    test.SetValley(this->Valley);
    test.SetManifold(this->Manifold);

    test.generateTriangles(RegionThreshold, AngleThreshold, MinDataValue, MaxDataValue, EigenValueThreshold, WaveSpeed, StencilRange);
	
    if (allocated_scalar_data) delete[] in_scalar;
    if (allocated_vector_data) delete[] in_grad;
    if (allocated_tensor_data) delete[] in_hesse;
   

    vtkSmartPointer<vtkPolyData> opd = test.GetPolyDataPointer();
    vtkSmartPointer<vtkImageData> oid = test.GetEigenValueGrid();

    vtkDataArray *dataArray = this->GetInputArrayToProcess(3, inputVector);
    if( dataArray )
    {
        int type = dataArray->GetDataType();
        vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
        arr->SetNumberOfComponents(dataArray->GetNumberOfComponents());
        arr->SetName(dataArray->GetName());

        // do some probing for additional data on ridges
        vtkSmartPointer<vtkPoints> ridgePoints = opd->GetPoints();

        for( vtkIdType ptId = 0; ptId < ridgePoints->GetNumberOfPoints(); ptId++ ) 
        {
            double pos[3], pcoords[3], weights[8];
            ridgePoints->GetPoint(ptId, pos);
            int subId=0;
            vtkCell *cell = in_grid->FindAndGetCell(pos, NULL, -1, 1e-5, subId, pcoords, weights);
        
            double data[9], result[9]={0,0,0,0,0,0,0,0,0};
            
            if(cell ) {
                for( vtkIdType cellPtId = 0; cellPtId<cell->GetNumberOfPoints(); cellPtId++)
                {
                    vtkIdType pt = cell->GetPointId(cellPtId);
                    dataArray->GetTuple(pt, data);

                    for( int i=0; i<dataArray->GetNumberOfComponents(); i++){
                        result[i] += data[i]*weights[cellPtId];
                    }
                }
            }
            arr->InsertNextTuple(result);
        }

        opd->GetPointData()->AddArray(arr);
        arr->Delete();
    }

    //send to pipeline
    output->ShallowCopy(opd);
    output1->ShallowCopy(oid);
    return 1;
}

void vtkRidge::getData(vtkSmartPointer<vtkImageData> in_grid, double* &in_data, int numComponents, bool &allocated_in_data)
{
    switch (in_grid->GetScalarType()){
    case VTK_FLOAT:
        {
            const size_t size = dim[0]*dim[1]*dim[2]*numComponents;
            std::vector<float> tmp_data(size);
            memcpy(tmp_data.data(), in_grid->GetScalarPointer(), size*sizeof(float));
            std::vector<double> tmp_double(tmp_data.begin(),tmp_data.end());
            in_data = new double[size];
            allocated_in_data = true;
            memcpy(in_data, tmp_double.data(), size*sizeof(double));
        }
        break;
    case VTK_DOUBLE:
    default:
        in_data = static_cast<double *>(in_grid->GetScalarPointer());
        break;
    }
}

////////// External Operators /////////////

void vtkRidge::PrintSelf(ostream &os, vtkIndent indent)
{
}
