#include "ioutils.h"
#include "morphutils.h"
#include "rjbutilities.h"

#include "itkMaskAdaptorImageFilter.h"
#include "itkTriangleThresholdImageFilter.h"
#include <itkOtsuThresholdImageFilter.h>
#include <itkMaskImageFilter.h>

#include <itkSmartPointer.h>
namespace itk
{
    template <typename T>
    class Instance : public T::Pointer {
    public:
        Instance() : SmartPointer<T>( T::New() ) {}
    };
}




int main(int argc, char * argv[])
{
  const unsigned dim = 3;
  typedef itk::Image<unsigned char, dim> LabImType;
  typedef itk::Image<float, dim> RawImType;

  RawImType::Pointer raw = readIm<RawImType>(argv[1]);

  itk::Instance <itk::OtsuThresholdImageFilter<RawImType, LabImType > > Otsu;
  Otsu->SetInput(raw);
  Otsu->SetOutsideValue(1);
  Otsu->SetInsideValue(0);

  writeIm<LabImType>(Otsu->GetOutput(), argv[2]);
  std::cout << "Otsu threshold: " << (float)Otsu->GetThreshold() << std::endl;
  
  itk::Instance <itk::TriangleThresholdImageFilter<RawImType, LabImType > > Tri;
  Tri->SetInput(raw);
  Tri->SetOutsideValue(1);
  Tri->SetInsideValue(0);

  writeIm<LabImType>(Tri->GetOutput(), argv[3]);
  std::cout << "Triangle threshold: " << (float)Tri->GetThreshold() << std::endl;
  
  return(EXIT_SUCCESS);
}

