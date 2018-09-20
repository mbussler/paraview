
#ifndef __GenerateBonds_h
#define __GenerateBonds_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"

class GenerateBonds : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(GenerateBonds,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct with user-specified implicit function; InsideOut turned off;
  // value set to 0.0; and generate clip scalars turned off.
  static GenerateBonds *New();

  vtkSetMacro(NeighborhoodSize,double);
  vtkGetMacro(NeighborhoodSize,double);

protected:
  GenerateBonds();
  ~GenerateBonds();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  double NeighborhoodSize;

private:

  GenerateBonds(const GenerateBonds&);  // Not implemented.
  void operator=(const GenerateBonds&);  // Not implemented.
};

#endif
