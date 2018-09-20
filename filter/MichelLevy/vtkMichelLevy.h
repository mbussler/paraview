/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMichelLevy.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __MichelLevy_h
#define __MichelLevy_h

#include "vtkDataSetAlgorithm.h"
#include "linalg.h"

class vtkMichelLevy : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkMichelLevy, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkMichelLevy *New();

    // Angle between the polarizer and the analyzer
    vtkSetClampMacro( AnalyzerAngle,double,0.0,180.0);
    vtkGetMacro( AnalyzerAngle,double);

    vtkSetClampMacro( ObjectAngle,double,-90.0,90.0);
    vtkGetMacro( ObjectAngle,double);

    vtkSetClampMacro( Thickness,double,0.0,100.0);
    vtkGetMacro( Thickness,double);

    vtkSetClampMacro( FringeValue,double,0.0,100.0);
    vtkGetMacro( FringeValue,double);

    vtkSetClampMacro( Wavelength,double,380.0,780.0);
    vtkGetMacro( Wavelength,double);

    vtkSetClampMacro( Intensity,double,0.001,10.0);
    vtkGetMacro( Intensity,double);

    vtkSetMacro( Monochromatic, bool);
    vtkGetMacro( Monochromatic, bool);
    
    vtkSetMacro( Lightfield, bool);
    vtkGetMacro( Lightfield, bool);

    vtkSetClampMacro( Gamma,double,0.001,10.0);
    vtkGetMacro( Gamma,double);
    
  protected:
    vtkMichelLevy();
    ~vtkMichelLevy();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkMichelLevy(const vtkMichelLevy&);  // Not implemented.
    void operator=(const vtkMichelLevy&);  // Not implemented.
  
    inline double TransFixed( double Gamma, double lambda);
    inline double Trans( double Gamma, double lambda, double alpha);
    inline double T( double delta, double lambda );
    inline double TL( double delta, double lambda );
    
    void createBifringenceTable();
    
    void getRGB( const double& retardation, vec3& rgb );
    void getRGB( const double& retardation, const double& alpha, vec3& rgb );
    
    vec3* RGB_lin;
    mat3 adobe_rgb;

    
    double AnalyzerAngle;
    double ObjectAngle;
    double Thickness;
    double FringeValue;
    double Wavelength;
    double Intensity;
    double Gamma;
    bool Monochromatic;
    bool Lightfield;
    
};

#endif
