#ifndef __itkKittlerIllingworthThresholdImageFilter_txx
#define __itkKittlerIllingworthThresholdImageFilter_txx
#include "itkKittlerIllingworthThresholdImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkKittlerIllingworthThresholdImageCalculator.h"
#include "itkProgressAccumulator.h"

namespace itk {

template<class TInputImage, class TOutputImage>
KittlerIllingworthThresholdImageFilter<TInputImage, TOutputImage>
::KittlerIllingworthThresholdImageFilter()
{
  m_OutsideValue   = NumericTraits<OutputPixelType>::Zero;
  m_InsideValue    = NumericTraits<OutputPixelType>::max();
  m_Threshold      = NumericTraits<InputPixelType>::Zero;
  m_NumberOfHistogramBins = 128;
}

template<class TInputImage, class TOutputImage>
void
KittlerIllingworthThresholdImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  typename ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  // Compute the KittlerIllingworth Threshold for the input image
  typename KittlerIllingworthThresholdImageCalculator<TInputImage>::Pointer calculator =
    KittlerIllingworthThresholdImageCalculator<TInputImage>::New();
  calculator->SetImage (this->GetInput());
  calculator->SetNumberOfHistogramBins (m_NumberOfHistogramBins);
  calculator->Compute();
  m_Threshold = calculator->GetThreshold();

  typename BinaryThresholdImageFilter<TInputImage,TOutputImage>::Pointer threshold =
    BinaryThresholdImageFilter<TInputImage,TOutputImage>::New();

  progress->RegisterInternalFilter(threshold,.5f);
  threshold->GraftOutput (this->GetOutput());
  threshold->SetInput (this->GetInput());
  threshold->SetLowerThreshold(NumericTraits<InputPixelType>::NonpositiveMin());
  threshold->SetUpperThreshold(calculator->GetThreshold());
  threshold->SetInsideValue (m_InsideValue);
  threshold->SetOutsideValue (m_OutsideValue);
  threshold->Update();

  this->GraftOutput(threshold->GetOutput());
}

template<class TInputImage, class TOutputImage>
void
KittlerIllingworthThresholdImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  TInputImage * input = const_cast<TInputImage *>(this->GetInput());
  if( input )
    {
    input->SetRequestedRegionToLargestPossibleRegion();
    }
}

template<class TInputImage, class TOutputImage>
void 
KittlerIllingworthThresholdImageFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "OutsideValue: "
     << static_cast<typename NumericTraits<OutputPixelType>::PrintType>(m_OutsideValue) << std::endl;
  os << indent << "InsideValue: "
     << static_cast<typename NumericTraits<OutputPixelType>::PrintType>(m_InsideValue) << std::endl;
  os << indent << "NumberOfHistogramBins: "
     << m_NumberOfHistogramBins << std::endl;
  os << indent << "Threshold (computed): "
     << static_cast<typename NumericTraits<InputPixelType>::PrintType>(m_Threshold) << std::endl;

}


}// end namespace itk
#endif
