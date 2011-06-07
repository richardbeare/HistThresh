#include "ioutils.h"
#include "itkTriangleThresholdImageFilter.h"

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
  
  itk::Instance <itk::TriangleThresholdImageFilter<RawImType, LabImType > > Tri;
  Tri->SetInput(raw);
  Tri->SetOutsideValue(1);
  Tri->SetInsideValue(0);
  Tri->DebugOn();
  writeIm<LabImType>(Tri->GetOutput(), argv[2]);
  std::cout << "Triangle threshold: " << (float)Tri->GetThreshold() << std::endl;
  
  return(EXIT_SUCCESS);
}

