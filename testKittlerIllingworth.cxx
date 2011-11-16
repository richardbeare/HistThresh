#include "ioutils.h"

#include "itkKittlerIllingworthThresholdImageFilter.h"

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

  
  itk::Instance <itk::KittlerIllingworthThresholdImageFilter<RawImType, LabImType > > Thr;
  Thr->SetInput(raw);
  Thr->SetOutsideValue(255);
  Thr->SetInsideValue(0);

  writeIm<LabImType>(Thr->GetOutput(), argv[2]);
  std::cout << "KittlerIllingworth threshold: " << (float)Thr->GetThreshold() << std::endl;
  

  return(EXIT_SUCCESS);
}

