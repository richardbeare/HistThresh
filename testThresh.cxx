#include "ioutils.h"
#include "morphutils.h"
#include "rjbutilities.h"

#include "itkMaskAdaptorImageFilter.h"
#include "itkTriangleThresholdImageFilter.h"
#include <itkOtsuThresholdImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkMaskNegatedImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>

#include <itkSmartPointer.h>
namespace itk
{
    template <typename T>
    class Instance : public T::Pointer {
    public:
        Instance() : SmartPointer<T>( T::New() ) {}
    };
}


template <class RawImType, class MaskImType>
typename RawImType::Pointer doNewMorphGrad(typename RawImType::Pointer raw, typename MaskImType::Pointer roughseg, double sigma2)
{
  itk::Instance< itk::MaskImageFilter <RawImType, MaskImType, RawImType> > Masker;
  Masker->SetInput(raw);
  Masker->SetInput2(roughseg);

  // 1 voxel dilation
  typename RawImType::Pointer g1 = doGradientOuter<RawImType>(Masker->GetOutput(), 1); 
  // mask again
  itk::Instance< itk::MaskImageFilter <RawImType, MaskImType, RawImType> > Masker2;
  Masker2->SetInput(g1);
  Masker2->SetInput2(roughseg);

  if (sigma2 > 0)
    {
    itk::Instance< itk::SmoothingRecursiveGaussianImageFilter< RawImType, RawImType > > Smoother;
    Smoother->SetInput(Masker2->GetOutput());
    Smoother->SetSigma(sigma2);
    typename RawImType::Pointer result = Smoother->GetOutput();
    result->Update();
    result->DisconnectPipeline();
    return(result);
    }
  else
    {
    typename RawImType::Pointer result = Masker2->GetOutput();
    result->Update();
    result->DisconnectPipeline();
    return(result);
    }

}


template <class LabImType, class RawImType>
void doRefine(typename LabImType::Pointer rough, typename LabImType::Pointer phase1Markers, typename RawImType::Pointer raw, int label, bool bright)
{
  // clipping of dark bits at the edge. Also using a watershed. Use
  // the first stage segmentation as a guide. Carry out an otsu
  // threshold on the stage 1 brain region and mask the dark section
  // out of the mask. Compute a gradient based control image and apply
  // watershed.

  typedef typename itk::MaskAdaptorImageFilter<RawImType, LabImType, LabImType> MaskAdaptorType;
  typename MaskAdaptorType::Pointer masker = MaskAdaptorType::New();

  typedef typename LabImType::PixelType PixelType;

  itk::Instance <itk::OtsuThresholdImageFilter< typename MaskAdaptorType::InternalInputImageType, typename MaskAdaptorType::InternalOutputImageType> > Otsu;
  Otsu->SetOutsideValue((PixelType)(bright));
  Otsu->SetInsideValue((PixelType)(!bright));
  masker->SetFilter(Otsu);
  masker->SetInput(raw);
  masker->SetMaskImage(rough);
  masker->SetDefaultValue(0);
  
  // We need to be very careful computing a gradient - we want to
  // select grey/csf transitions.

  // step 1: get the raw gradient
  typename RawImType::Pointer grad = doNewMorphGrad<RawImType, LabImType>(raw, rough, 0.0);

  // step 2: threshold the gradient within the rough mask
  // ranks seem the appropriate way to select the threshold
  std::vector<float> quantiles, qvals;
  quantiles.push_back(0.5);
  qvals = computeMaskQuantiles<RawImType, LabImType, float>(grad, rough, quantiles);
  std::cout << qvals[0] << std::endl;
  itk::Instance< itk::BinaryThresholdImageFilter<RawImType, LabImType> > RThreshType;
  RThreshType->SetInput(grad);
  RThreshType->SetUpperThreshold(qvals[0]);
  RThreshType->SetInsideValue(0);
  RThreshType->SetOutsideValue(1);


  typedef typename LabImType::PixelType PixelType;

  writeIm<LabImType>(RThreshType->GetOutput(), "/tmp/gradmask.nii.gz");
  writeIm<RawImType>(grad, "/tmp/grad.nii.gz");

  // now generate a mask of dark areas inside the edges
  typename MaskAdaptorType::Pointer maskerGrad = MaskAdaptorType::New();
  itk::Instance <itk::TriangleThresholdImageFilter< typename MaskAdaptorType::InternalInputImageType, typename MaskAdaptorType::InternalOutputImageType> > TriGrad;
  TriGrad->SetOutsideValue((PixelType)(bright));
  TriGrad->SetInsideValue((PixelType)(!bright));
  maskerGrad->SetFilter(TriGrad);
  maskerGrad->SetInput(raw);
  maskerGrad->SetMaskImage(RThreshType->GetOutput());
  maskerGrad->SetDefaultValue(0);
  maskerGrad->Update();
  std::cout << " threshold value " << (int)TriGrad->GetThreshold() << std::endl;
}
int main(int argc, char * argv[])
{
  const unsigned dim = 3;
  typedef itk::Image<unsigned char, dim> LabImType;
  typedef itk::Image<float, dim> RawImType;

  LabImType::Pointer rough = readIm<LabImType>("/tmp/rough.nii.gz");
  LabImType::Pointer segPhase1 = readIm<LabImType>("/tmp/segPhase1.nii.gz");
  RawImType::Pointer raw = readIm<RawImType>("/tmp/raw.nii.gz");

  doRefine<LabImType, RawImType>(rough, segPhase1, raw, 1, false);

}

